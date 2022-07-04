#pragma once

namespace hexa
{
	enum class InterpolationType { Drop, Linear, Lagrange3, BSpline3, CatmullRom, Opti3, Opti4 };

	template <typename Type, InterpolationType type>
	struct Interpolator;

	template <typename Type>
	struct Interpolator<Type, InterpolationType::Drop>
	{
		//==============================================================================
		static constexpr size_t order = 0;

		static constexpr size_t numPoints = 1;

		//==============================================================================
		Type operator() ([[maybe_unused]] double x, Type y0) const noexcept
		{
			return y0;
		}
	};

	template <typename Type>
	struct Interpolator<Type, InterpolationType::Linear>
	{
		//==============================================================================
		static constexpr size_t order = 1;

		static constexpr size_t numPoints = 2;

		//==============================================================================
		Type operator() (double x, Type y0, Type y1) const noexcept
		{
			return static_cast<Type>(y1 + x * (y0 - y1));
		}
	};

	template <typename Type>
	struct Interpolator<Type, InterpolationType::BSpline3>
	{
		//==============================================================================
		static constexpr size_t order = 3;

		static constexpr size_t numPoints = 4;

		//==============================================================================
		Type operator()(Type x, Type y0, Type y1, Type y2, Type y3) const noexcept
		{
			const double a = 0.1666666667 * y0 + 0.6666666667 * y1 + 0.1666666667 * y2,
				b = -0.5 * y0 + 0.5 * y2,
				c = 0.5 * y0 - y1 + 0.5 * y2,
				d = -0.1666666667 * y0 + 0.5 * y1 - 0.5 * y2 + 0.1666666667 * y3;
			return static_cast<Type>(a + x * (b + x * (c + d * x)));;
		}
	};

	template <typename Type>
	struct Interpolator<Type, InterpolationType::CatmullRom>
	{
		//==============================================================================
		static constexpr size_t order = 3;

		static constexpr size_t numPoints = 4;

		//==============================================================================
		Type operator() (double x, Type y0, Type y1, Type y2, Type y3) const noexcept
		{
			const double a = y1,
				b = -0.5 * y0 + 0.5 * y2,
				c = y0 - 2.5 * y1 + 2. * y2 - 0.5 * y3,
				d = -0.5 * y0 + 1.5 * y1 - 1.5 * y2 + 0.5 * y3;
			return static_cast<Type>(a + x * (b + x * (c + d * x)));
		}
	};

	template <typename Type>
	struct Interpolator<Type, InterpolationType::Lagrange3>
	{
		//==============================================================================
		static constexpr size_t order = 3;

		static constexpr size_t numPoints = 4;

		//==============================================================================
		Type operator() (double x, Type y0, Type y1, Type y2, Type y3) const noexcept
		{
			const double a = y1,
				b = -0.3333333333 * y0 - 0.5 * y1 + y2 - 0.1666666667 * y3,
				c = 0.5 * y0 - y1 + 0.5 * y2,
				d = -0.1666666667 * y0 + 0.5 * y1 - 0.5 * y2 + 0.1666666667 * y3;
			return static_cast<Type>(a + x * (b + x * (c + d * x)));
		}
	};

	template <typename Type>
	struct Interpolator<Type, InterpolationType::Opti3>
	{
		//==============================================================================
		static constexpr size_t order = 3;

		static constexpr size_t numPoints = 4;

		//==============================================================================
		Type operator() (double x, Type y0, Type y1, Type y2, Type y3) const noexcept
		{
			const double a = 0.20345744715566433 * y0 + 0.5924449242027232 * y1 + 0.20184198969656253 * y2 + 0.002240727070748738 * y3,
				b = -0.49823192036183106 * y0 + 0.03573669883299365 * y1 + 0.45663331520682054 * y2 + 0.005951377567825489 * y3,
				c = 0.3987650580367404 * y0 - 0.7866488859776489 * y1 + 0.29427887193783475 * y2 + 0.09351548475726523 * y3,
				d = -0.10174985775982505 * y0 + 0.36030925263849456 * y1 - 0.36030925263849456 * y2 + 0.10174985775982505 * y3;
			return static_cast<Type>(a + x * (b + x * (c + d * x)));
		}
	};

	template <typename Type>
	struct Interpolator<Type, InterpolationType::Opti4>
	{
		//==============================================================================
		static constexpr size_t order = 4;

		static constexpr size_t numPoints = 4;

		//==============================================================================
		Type operator() (double x, Type y0, Type y1, Type y2, Type y3) const noexcept
		{
			const double a = 0.1882882705853767 * y0 + 0.6234244946593812 * y1 + 0.18828618868306463 * y2 + 1.0459580992439044e-6 * y3,
				b = -0.521495184193515 * y0 + 0.06456923251842608 * y1 + 0.43534733489266775 * y2 + 0.021578619812177235 * y3,
				c = 0.4995463562127711 * y0 - 0.9990450958317605 * y1 + 0.49945043290341395 * y2 + 0.00004829195935707187 * y3,
				d = -0.16617743854192843 * y0 + 0.49917660509564427 * y1 - 0.49982041406113864 * y2 + 0.16682127096034755 * y3,
				e = -0.00016095810460478 * y0 + 0.0001609522413736 * y1 + 0.0001609522413736 * y2 - 0.00016095810460478 * y3;
			return static_cast<Type>(a + x * (b + x * (c + x * (d + e * x))));
		}
	};
}