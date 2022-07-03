#pragma once

namespace hexa
{
	template <typename T>
	struct c
	{
		/** A predefined value for Pi */
		static constexpr T pi
			= static_cast<T>(3.141592653589793238462643383279502884197169399375105820974944592307816L);

		/** A predefined value for 2 * Pi */
		static constexpr T twoPi
			= static_cast<T>(6.283185307179586476925286766559005768394338798750211641949889184615633L);

		/** A predefined value for 2 * Pi */
		static constexpr T reciprPi
			= static_cast<T>(0.318309886183790671537767526745028724068919291480912897495334688117793L);

		/** A predefined value for 1 / (2 * Pi) */
		static constexpr T reciprTwoPi
			= static_cast<T>(0.159154943091895335768883763372514362034459645740456448747667344058896L);

		/** A predefined value for Pi / 2 */
		static constexpr T halfPi
			= static_cast<T>(1.570796326794896619231321691639751442098584699687552910487472296153908L);

		/** A predefined value for Euler's number */
		static constexpr T e
			= static_cast<T>(2.718281828459045235360287471352662497757247093699959574966967627724077L);

		/** A predefined value for square of Euler's number */
		static constexpr T eSquare
			= static_cast<T>(7.389056098930650227230427460575007813180315570551847324087127822522574L);

		/** A predefined value for square of Pi number */
		static constexpr T piSquare
			= static_cast<T>(9.869604401089358618834490999876151135313699407240790626413349376220045L);

		/** A predefined value for sqrt(2) */
		static constexpr T sqrt2
			= static_cast<T>(1.414213562373095048801688724209698078569671875376948073176679737990732L);

		/** A predefined value for 1/sqrt(2) */
		static constexpr T reciprSqrt2
			= static_cast<T>(0.7071067811865475244008443621048490392848359376884740365883398689953662L);

		/** A predefined value for square root of 3 */
		static constexpr T sqrt3
			= static_cast<T>(1.732050807568877293527446341505872366942805253810380628055806979451933L);

		/** A predefined value for square root of 5 */
		static constexpr T sqrt5
			= static_cast<T>(2.236067977499789696409173668731276235440618359611525724270897245410521L);

		/** A predefined value for the 2 in 3/4 power */
		static constexpr T two35
			= static_cast<T>(1.681792830507429086062250952466429790080068524713569021626452171949850L);

		/** A predefined value for Pi/4 */
		static constexpr T quarterPi
			= static_cast<T>(0.7853981633974483096156608458198757210492923498437764552437361480769541L);

		/** A predefined value for Pi^4 */
		static constexpr T piPow4
			= static_cast<T>(97.40909103400243723644033268870511124972758567268542169146785938997086L);

		/** A predefined value for Ln(2) */
		static constexpr T ln2
			= static_cast<T>(0.6931471805599453094172321214581765680755001343602552541206800094933936L);

		/** A predefined value for Ln(3) */
		static constexpr T ln3
			= static_cast<T>(1.098612288668109691395245236922525704647490557822749451734694333637494L);

		/** A predefined value for Ln(6) */
		static constexpr T ln6
			= static_cast<T>(1.791759469228055000812477358380702272722990692183004705855374343130888L);

		/** A predefined value for Ln(10) */
		static constexpr T ln10
			= static_cast<T>(2.302585092994045684017991454684364207601101488628772976033327900967573L);

		/** A predefined value for Euler–Mascheroni constant */
		static constexpr T eulerGamma
			= static_cast<T>(0.5772156649015328606065120900824024310421593359399235988057672348848677L);

		/** A predefined value for golden ratio constant */
		static constexpr T goldenRatio
			= static_cast<T>(1.61803398874989484820458683436563811772030917980576286213544862270526L);

		/** A predefined value for silver ratio constant */
		static constexpr T silverRatio
			= static_cast<T>(2.414213562373095048801688724209698078569671875376948073176679737990732L);

		/** A predefined value for bronze ratio constant */
		static constexpr T bronzeRatio
			= static_cast<T>(3.302775637731994646559610633735247973125648286922623106355226528113583L);

		/** A predefined value for Catalan's constant */
		static constexpr T catalan
			= static_cast<T>(0.915965594177219015054603514932384110774149374281672134266498119621763L);

		/** A predefined value for Apery's constant */
		static constexpr T zeta3
			= static_cast<T>(1.202056903159594285399738161511449990764986292340498881792271555341838L);

		/** A predefined value for Khinchin's constant */
		static constexpr T khc
			= static_cast<T>(2.685452001065306445309714835481795693820382293994462953051152345557219L);

		/** A predefined value for Ramanujan–Soldner constant */
		static constexpr T mu
			= static_cast<T>(1.451369234883381050283968485892027449493032283648015863093004557662426L);

		/** A predefined value for log2(e) */
		static constexpr T log2e
			= static_cast<T>(1.442695040888963407359924681001892137426645954152985934135449406931109L);

		/** A predefined value for log2(10) */
		static constexpr T log2p10
			= static_cast<T>(3.321928094887362347870319429489390175864831393024580612054756395815934L);

		/** A predefined value for reciprocal of log2(e) */
		static constexpr T reciprLog2e
			= static_cast<T>(0.693147180559945309417232121458176568075500134360255254120680009493394L);

		/** A predefined value for reciprocal of log2(10) */
		static constexpr T reciprLog2p10
			= static_cast<T>(0.3010299956639811952137388947244930267681898814621085413104274611271083L);
	};

	using cd = c<double>;
	using cf = c<float>;
}