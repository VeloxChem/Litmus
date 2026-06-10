// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.
// E-mail: rinkevic@kth.se
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "config.hpp"

#include <cctype>
#include <fstream>
#include <sstream>

namespace cfg {  // cfg namespace

ConfigError::ConfigError(const std::string& message)

    : std::runtime_error(message)
{
}

void
Config::set(const std::string& key, Value value)
{
    _values[key] = std::move(value);
}

bool
Config::has(const std::string& key) const
{
    return _values.find(key) != _values.end();
}

const Config::Value&
Config::_at(const std::string& key) const
{
    const auto it = _values.find(key);

    if (it == _values.end())
    {
        throw ConfigError("config: missing required key '" + key + "'");
    }

    return it->second;
}

std::string
Config::get_string(const std::string& key) const
{
    const auto& value = _at(key);

    if (const auto* str = std::get_if<std::string>(&value))
    {
        return *str;
    }

    throw ConfigError("config: key '" + key + "' is not a string");
}

std::string
Config::get_string(const std::string& key, const std::string& fallback) const
{
    return has(key) ? get_string(key) : fallback;
}

int
Config::get_int(const std::string& key) const
{
    const auto& value = _at(key);

    if (const auto* num = std::get_if<int>(&value))
    {
        return *num;
    }

    throw ConfigError("config: key '" + key + "' is not an integer");
}

int
Config::get_int(const std::string& key, int fallback) const
{
    return has(key) ? get_int(key) : fallback;
}

bool
Config::get_bool(const std::string& key) const
{
    const auto& value = _at(key);

    if (const auto* flag = std::get_if<bool>(&value))
    {
        return *flag;
    }

    throw ConfigError("config: key '" + key + "' is not a boolean");
}

bool
Config::get_bool(const std::string& key, bool fallback) const
{
    return has(key) ? get_bool(key) : fallback;
}

std::vector<int>
Config::get_int_array(const std::string& key) const
{
    const auto& value = _at(key);

    if (const auto* arr = std::get_if<std::vector<int>>(&value))
    {
        return *arr;
    }

    throw ConfigError("config: key '" + key + "' is not an integer array");
}

std::vector<int>
Config::get_int_array(const std::string& key, const std::vector<int>& fallback) const
{
    return has(key) ? get_int_array(key) : fallback;
}

namespace {  // unnamed namespace for parsing helpers

/// Removes leading and trailing ASCII whitespace.
std::string
trim(const std::string& text)
{
    std::size_t beg = 0;

    std::size_t end = text.size();

    while ((beg < end) && std::isspace(static_cast<unsigned char>(text[beg]))) beg++;

    while ((end > beg) && std::isspace(static_cast<unsigned char>(text[end - 1]))) end--;

    return text.substr(beg, end - beg);
}

/// Drops a trailing '#' comment, honouring '#' inside a double-quoted string.
std::string
strip_comment(const std::string& line)
{
    bool in_quotes = false;

    for (std::size_t i = 0; i < line.size(); i++)
    {
        const auto ch = line[i];

        if (ch == '"') in_quotes = !in_quotes;

        if ((ch == '#') && !in_quotes) return line.substr(0, i);
    }

    return line;
}

/// Checks whether a token is a signed decimal integer.
bool
is_integer(const std::string& token)
{
    if (token.empty()) return false;

    std::size_t pos = ((token[0] == '+') || (token[0] == '-')) ? 1 : 0;

    if (pos == token.size()) return false;

    for (std::size_t i = pos; i < token.size(); i++)
    {
        if (!std::isdigit(static_cast<unsigned char>(token[i]))) return false;
    }

    return true;
}

/// Parses a single integer token at the given line number.
int
parse_int(const std::string& token, std::size_t lineno)
{
    if (!is_integer(token))
    {
        throw ConfigError("config: line " + std::to_string(lineno) +
                          ": expected an integer, got '" + token + "'");
    }

    return std::stoi(token);
}

/// Parses a bracketed, comma-separated integer array.
std::vector<int>
parse_array(const std::string& raw, std::size_t lineno)
{
    if (raw.back() != ']')
    {
        throw ConfigError("config: line " + std::to_string(lineno) +
                          ": unterminated array '" + raw + "'");
    }

    const auto inner = trim(raw.substr(1, raw.size() - 2));

    std::vector<int> values;

    if (inner.empty()) return values;

    std::stringstream stream(inner);

    std::string token;

    while (std::getline(stream, token, ','))
    {
        values.push_back(parse_int(trim(token), lineno));
    }

    return values;
}

/// Parses the right-hand side of an assignment into a typed value.
Config::Value
parse_value(const std::string& raw, std::size_t lineno)
{
    if (raw.empty())
    {
        throw ConfigError("config: line " + std::to_string(lineno) + ": missing value");
    }

    if (raw.front() == '[') return parse_array(raw, lineno);

    if (raw.front() == '"')
    {
        if ((raw.size() < 2) || (raw.back() != '"'))
        {
            throw ConfigError("config: line " + std::to_string(lineno) +
                              ": unterminated string " + raw);
        }

        return raw.substr(1, raw.size() - 2);
    }

    if (raw == "true") return true;

    if (raw == "false") return false;

    if (is_integer(raw)) return std::stoi(raw);

    // anything else is taken as a bare (unquoted) string

    return raw;
}

}  // namespace

Config
parse_string(const std::string& text)
{
    Config config;

    std::stringstream stream(text);

    std::string line;

    std::size_t lineno = 0;

    while (std::getline(stream, line))
    {
        lineno++;

        const auto content = trim(strip_comment(line));

        if (content.empty()) continue;

        const auto eq = content.find('=');

        if (eq == std::string::npos)
        {
            throw ConfigError("config: line " + std::to_string(lineno) +
                              ": expected 'key = value', got '" + content + "'");
        }

        const auto key = trim(content.substr(0, eq));

        if (key.empty())
        {
            throw ConfigError("config: line " + std::to_string(lineno) + ": empty key");
        }

        config.set(key, parse_value(trim(content.substr(eq + 1)), lineno));
    }

    return config;
}

Config
parse_file(const std::string& path)
{
    std::ifstream stream(path);

    if (!stream.is_open())
    {
        throw ConfigError("config: cannot open file '" + path + "'");
    }

    std::stringstream buffer;

    buffer << stream.rdbuf();

    return parse_string(buffer.str());
}

}  // namespace cfg
