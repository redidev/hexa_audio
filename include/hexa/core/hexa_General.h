#pragma once

#include <limits>
#include <type_traits>

namespace hexa::utils
{
	//==============================================================================
/** Finds lowest near value, that is a GCD of a factor (factor needs to be a power of 2). */
	template <typename T, typename = std::enable_if_t<std::is_integral_v<std::remove_cv_t<T>>>>
	constexpr T roundDownToMultiple(T value, T factor) { return value & -factor; }

	/** Finds next power of 2, greater than value. */
	template <typename T, typename = std::enable_if_t<std::is_integral_v<std::remove_cv_t<T>>>>
	constexpr T nextPowerOfTwo(T value)
	{
		--value;
		value |= value >> 1;
		value |= value >> 2;
		value |= value >> 4;
		value |= value >> 8;
		value |= value >> 16; 
		++value;
		return value;
	}

	//==============================================================================
	/** Finds lowest near x, that is a multiple of a factor (factor needs to be a power of 2). (See Bill Prisoner's 'From green to red' book) */
	template <typename Type, typename = std::enable_if_t<std::is_integral_v<std::remove_cv_t<Type>>>>
	constexpr Type alignDown(Type x, Type align) { return x & ~(align - 1); }

	/** Finds greatest near value, that is a multiple of a factor (factor needs to be a power of 2). (See Bill Prisoner's book) */
	template <typename Type, typename = std::enable_if_t<std::is_integral_v<std::remove_cv_t<Type>>>>
	constexpr Type alignUp(Type x, Type align)
	{
		const auto alignMinus1 = align - 1;
		return (x + alignMinus1) & (~alignMinus1);
	}

	//==============================================================================
	/** Remaps a value from a 0..1 range to a target range. */
	template <typename Type, typename = std::enable_if_t<std::is_floating_point_v<std::remove_cv_t<Type>>>>
	constexpr Type map(Type value0to1, Type leftValue, Type rightValue)
	{
		return leftValue + (rightValue - leftValue) * value0to1;
	}

	/** Remaps a value from a 0..1 range to a target range 0..rightValue. */
	template <typename Type, typename = std::enable_if_t<std::is_floating_point_v<std::remove_cv_t<Type>>>>
	constexpr Type map(Type value0to1, Type rightValue)
	{
		return rightValue * value0to1;
	}

	/** Remaps a value from a source range to a target range. */
	template <typename Type, typename = std::enable_if_t<std::is_floating_point_v<std::remove_cv_t<Type>>>>
	constexpr Type map(Type valueToMap, Type leftValueSrc, Type rightValueSrc, Type leftValueDst, Type rightValueDst)
	{
		return leftValueDst + ((rightValueDst - leftValueDst) * (valueToMap - leftValueSrc)) / (rightValueSrc - leftValueSrc);
	}

	//==============================================================================
	/** Checks two floating-point numbers on equality. */
	template <typename Type, typename = std::enable_if_t<std::is_floating_point_v<std::remove_cv_t<Type>>>>
	constexpr bool areSame(const Type& value1, const Type& value2)
	{
		return std::abs(value1 - value2) < std::numeric_limits<Type>::epsilon();
	}

	// FNV-1a hash, 32-bit 
	constexpr std::uint32_t fnv1a(const char* str, std::uint32_t hash = 2166136261UL)
	{
		return *str ? fnv1a(str + 1, (hash ^ *str) * 16777619ULL) : hash;
	}
}