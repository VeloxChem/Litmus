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

#ifndef config_hpp
#define config_hpp

#include <map>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace cfg {  // cfg namespace

/// Error thrown on malformed configuration input or invalid value access.
class ConfigError : public std::runtime_error
{
public:
    explicit ConfigError(const std::string& message);
};

/// A parsed run configuration: a flat map of keys to typed values.
///
/// Supported value types are string, integer, boolean, and integer array.
/// Values are read back through the typed accessors, each of which has an
/// overload taking a fallback used when the key is absent.
class Config
{
public:
    /// The value variant stored against each key.
    using Value = std::variant<std::string, int, bool, std::vector<int>>;

    /// Inserts or replaces the value bound to a key.
    /// @param key The configuration key.
    /// @param value The value to bind.
    void set(const std::string& key, Value value);

    /// Checks whether a key is present.
    /// @param key The configuration key.
    /// @return True if the key is present, False otherwise.
    bool has(const std::string& key) const;

    /// Reads a string value.
    /// @param key The configuration key.
    /// @return The string value (throws ConfigError if absent or not a string).
    std::string get_string(const std::string& key) const;

    /// Reads a string value, falling back when the key is absent.
    /// @param key The configuration key.
    /// @param fallback The value to return when the key is absent.
    /// @return The string value, or the fallback.
    std::string get_string(const std::string& key, const std::string& fallback) const;

    /// Reads an integer value.
    /// @param key The configuration key.
    /// @return The integer value (throws ConfigError if absent or not an integer).
    int get_int(const std::string& key) const;

    /// Reads an integer value, falling back when the key is absent.
    /// @param key The configuration key.
    /// @param fallback The value to return when the key is absent.
    /// @return The integer value, or the fallback.
    int get_int(const std::string& key, int fallback) const;

    /// Reads a boolean value.
    /// @param key The configuration key.
    /// @return The boolean value (throws ConfigError if absent or not a boolean).
    bool get_bool(const std::string& key) const;

    /// Reads a boolean value, falling back when the key is absent.
    /// @param key The configuration key.
    /// @param fallback The value to return when the key is absent.
    /// @return The boolean value, or the fallback.
    bool get_bool(const std::string& key, bool fallback) const;

    /// Reads an integer-array value.
    /// @param key The configuration key.
    /// @return The integer array (throws ConfigError if absent or not an array).
    std::vector<int> get_int_array(const std::string& key) const;

    /// Reads an integer-array value, falling back when the key is absent.
    /// @param key The configuration key.
    /// @param fallback The value to return when the key is absent.
    /// @return The integer array, or the fallback.
    std::vector<int> get_int_array(const std::string& key, const std::vector<int>& fallback) const;

private:
    /// Returns the stored value for a key, throwing ConfigError if absent.
    const Value& _at(const std::string& key) const;

    /// The parsed key/value bindings.
    std::map<std::string, Value> _values;
};

/// Parses configuration text (the contents of a config file).
/// @param text The configuration text.
/// @return The parsed configuration (throws ConfigError on a malformed line).
Config parse_string(const std::string& text);

/// Parses a configuration file.
/// @param path The path to the configuration file.
/// @return The parsed configuration (throws ConfigError if the file cannot be
///         opened or contains a malformed line).
Config parse_file(const std::string& path);

}  // namespace cfg

#endif /* config_hpp */
