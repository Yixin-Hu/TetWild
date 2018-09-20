// Reenable the warnings disabled before including third party code
#if defined(__clang__)
#pragma clang diagnostic pop // -Wshadow
#pragma clang diagnostic pop // -Wsign-compare
#pragma clang diagnostic pop // -Wswitch-default
#pragma clang diagnostic pop // -Wformat-nonliteral
#pragma clang diagnostic pop // -Wswitch-enum
#pragma clang diagnostic pop // -Wstrict-overflow
// #pragma clang diagnostic pop // -Wnoexcept
#pragma clang diagnostic pop // -Wctor-dtor-privacy
#pragma clang diagnostic pop // -Wnull-dereference
#pragma clang diagnostic pop // -Wcast-qual
#pragma clang diagnostic pop // -Wmissing-noreturn
#pragma clang diagnostic pop // -Woverloaded-virtual
#pragma clang diagnostic pop // "-Wsign-promo"
#pragma clang diagnostic pop // "-Wcast-align"
#pragma clang diagnostic pop // "-Wnull-pointer-arithmetic"
#pragma clang diagnostic pop // "-Wc++17-extensions"
#elif (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop // -Wshadow
#pragma GCC diagnostic pop // -Wsign-compare
#pragma GCC diagnostic pop // -Wswitch-default
#pragma GCC diagnostic pop // -Wformat-nonliteral
#pragma GCC diagnostic pop // -Wswitch-enum
#pragma GCC diagnostic pop // -Wstrict-overflow
#pragma GCC diagnostic pop // -Wnoexcept
#pragma GCC diagnostic pop // -Wctor-dtor-privacy
#pragma GCC diagnostic pop // -Wnull-dereference
#pragma GCC diagnostic pop // -Wcast-qual
#pragma GCC diagnostic pop // -Wmissing-noreturn
#pragma GCC diagnostic pop // -Woverloaded-virtual
#pragma GCC diagnostic pop // "-Wsign-promo"
#endif
