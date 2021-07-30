#pragma once

#include <concepts>
#include <stdexcept>
#include <string>

#include <seqan/sequence.h>

template <typename in_t>
std::string to_str(in_t && in)
{
    return std::to_string(std::forward<in_t>(in));
}

template <typename in_t>
std::string to_str(in_t && in) requires std::constructible_from<std::string, in_t>
{
    return std::string(std::forward<in_t>(in));
}

inline std::string to_str(char in)
{
    return std::string{in};
}

inline std::string to_str(std::string in)
{
    return in;
}

inline std::string to_str(seqan::CharString && in)
{
    return seqan::toCString(in);
}

struct error : public std::runtime_error
{
    template <typename... ts>
    error(ts... what) : std::runtime_error{(to_str(std::forward<ts>(what)) + ...)}
    {}
};
