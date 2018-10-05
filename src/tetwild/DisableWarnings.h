// Disable compiler warnings before including third party code
#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wshadow"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wsign-compare"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wswitch-default"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wformat-nonliteral"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wswitch-enum"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wstrict-overflow"
	// #pragma clang diagnostic push
	// #pragma clang diagnostic ignored "-Wnoexcept"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wctor-dtor-privacy"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wnull-dereference"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wcast-qual"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wmissing-noreturn"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Woverloaded-virtual"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wsign-promo"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wcast-align"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wnull-pointer-arithmetic"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wc++17-extensions"
#elif (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wshadow"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wsign-compare"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wswitch-default"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wformat-nonliteral"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wswitch-enum"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wstrict-overflow"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wnoexcept"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wnull-dereference"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wcast-qual"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wmissing-noreturn"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Woverloaded-virtual"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wsign-promo"
#endif
