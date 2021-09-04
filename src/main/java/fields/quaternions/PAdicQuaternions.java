package fields.quaternions;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.local.PAdicField;
import fields.local.PAdicField.PAdicNumber;
import fields.local.Value;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;
import util.MiscAlgorithms;

public class PAdicQuaternions extends AbstractQuaternions<PAdicNumber> {
	private PAdicField base;
	private Reals reals = Reals.r(1024);
	private int hilbertSymbol;

	private PAdicQuaternions(PAdicField base, PAdicNumber a, PAdicNumber b) {
		super(base, a, b);
		this.base = base;
		this.hilbertSymbol = 0;
	}

	public static PAdicQuaternions quaternions(PAdicField base) {
		if (base.reduction().characteristic().equals(BigInteger.TWO)) {
			return new PAdicQuaternions(base, base.getInteger(-3), base.uniformizer());
		}
		return new PAdicQuaternions(base, base.lift(base.reduction().quadraticNonResidue()), base.uniformizer());
	}

	public static PAdicQuaternions quaternions(PAdicField base, PAdicNumber a, PAdicNumber b) {
		return new PAdicQuaternions(base, a, b);
	}

	@Override
	public int hilbertSymbol() {
		if (hilbertSymbol == 0) {
			int aIndex = upToSquaresIndex(base, a());
			int bIndex = upToSquaresIndex(base, b());
			int[][] result;
			if (base.reduction().characteristic().equals(BigInteger.TWO)) {
				result = new int[][] { { 1, -1, -1, 1, 1, -1, -1, 1 }, // -6
						{ -1, 1, -1, 1, 1, -1, 1, -1 }, // -3
						{ -1, -1, -1, -1, 1, 1, 1, 1 }, // -2
						{ 1, 1, -1, -1, 1, 1, -1, -1 }, // -1
						{ 1, 1, 1, 1, 1, 1, 1, 1 }, // 1
						{ -1, -1, 1, 1, 1, 1, -1, -1 }, // 2
						{ -1, 1, 1, -1, 1, -1, -1, 1 }, // 3
						{ 1, -1, 1, -1, 1, -1, 1, -1 } // 6
				};
			} else {
				int s = MiscAlgorithms.jacobiSymbol(BigInteger.valueOf(-1), base.reduction().characteristic());
				result = new int[][] { { 1, 1, 1, 1 }, // 1
						{ 1, 1, -1, -1 }, // e
						{ 1, -1, s, -s }, // p
						{ 1, -1, -s, s } // ep
				};
			}
			hilbertSymbol = result[aIndex][bIndex];
		}
		return hilbertSymbol;
	}

	/**
	 * Index of t in [-6, -3, -2, -1, 1, 2, 3, 6] if char = (0, 2), index of t in
	 * [1, e, p, ep] if char = (0, p), e non square residue
	 */
	private static int upToSquaresIndex(PAdicField base, PAdicNumber t) {
		if (base.reduction().characteristic().equals(BigInteger.TWO)) {
			if (t.equals(base.zero())) {
				return 4;
			}
			int value = base.valuation(t).value();
			t = base.rightShift(t, value);
			int mod8 = t.digit(0).intValueExact() + 2 * t.digit(1).intValueExact() + 4 * t.digit(2).intValueExact();
			if (value % 2 == 0) {
				switch (mod8) {
				case 1:
					return 4;
				case 3:
					return 6;
				case 5:
					return 1;
				case 7:
					return 3;
				default:
					throw new ArithmeticException("Number not even!");
				}
			} else {
				switch (mod8) {
				case 1:
					return 5;
				case 3:
					return 7;
				case 5:
					return 0;
				case 7:
					return 2;
				default:
					throw new ArithmeticException("Number not even!");
				}

			}
		}
		if (t.equals(base.zero())) {
			return 0;
		}
		int value = base.valuation(t).value();
		t = base.multiply(t, base.power(base.uniformizer(), -value));
		int sqrt = base.reduction().hasSqrt(base.reduce(t)) ? 0 : 1;
		value = value > 0 ? (value % 2) : ((-value) % 2);
		return value * 2 + sqrt;
	}

	@Override
	public PAdicNumber discriminant() {
		if (isIntegral()) {
			return base.uniformizer();
		}
		return base.one();
	}

	@Override
	public boolean isIntegral() {
		return hilbertSymbol() == -1;
	}

	public Value valuation(Quaternion<PAdicNumber> t) {
		return base.valuation(reducedNorm(t));
	}

	public Real value(Quaternion<PAdicNumber> t) {
		return reals.positiveSqrt(base.value(reducedNorm(t)));
	}

	private Value goodSplittingElement(Quaternion<PAdicNumber> candidate) {
		Quaternion<PAdicNumber> second = negativeCommutor(candidate);
		Quaternion<PAdicNumber> third = multiply(candidate, second);
		List<Vector<PAdicNumber>> reducedVectors = new ArrayList<>();
		reducedVectors.add(candidate.imaginaryPart());// reduceImaginaryVector(candidate));
		reducedVectors.add(second.imaginaryPart());// reduceImaginaryVector(second));
		reducedVectors.add(third.imaginaryPart());// reduceImaginaryVector(third));
		MatrixAlgebra<PAdicNumber> m = new FreeModule<>(base, 3).matrixAlgebra();
		return base.valuation(m.determinant(Matrix.fromColumns(reducedVectors)));
	}

	@Override
	protected Quaternion<PAdicNumber> normalize(Quaternion<PAdicNumber> t) {
		if (t.equals(zero())) {
			return zero();
		}
		PAdicNumber nonZero = null;
		Value minValue = Value.INFINITY;
		Vector<PAdicNumber> asVector = asVector(t);
		boolean allEqual = true;
		for (int i = 0; i < 4; i++) {
			Value value = base.valuation(asVector.get(i + 1));
			if (value.isInfinite()) {
				continue;
			}
			if (nonZero == null) {
				nonZero = asVector.get(i + 1);
			} else 
			if (allEqual && !asVector.get(i+1).equals(nonZero) && !asVector.get(i+1).equals(base.negative(nonZero))) {
				allEqual = false;
			}
			if (minValue.compareTo(value) > 0) {
				minValue = value;
			}
		}
		PAdicNumber divisor = allEqual ? nonZero : base.one(); // base.leftShift(base.one(), minValue.value());
		return divide(t, getEmbedding(divisor));
	}

	@Override
	protected Quaternion<PAdicNumber> splittingElement() {
		PrimeField reduction = base.reduction();
		int aIndex = upToSquaresIndex(base, a());
		int aValueSqrt = MiscAlgorithms.DivRoundDown(base.valuation(a()).value(), 2);
		int bIndex = upToSquaresIndex(base, b());
		int bValueSqrt = MiscAlgorithms.DivRoundDown(base.valuation(b()).value(), 2);
		PAdicNumber x;
		PAdicNumber y;
		if (reduction.characteristic().equals(BigInteger.TWO)) {
			if (bIndex == 4) {
				x = base.one();
				y = base.zero();
			} else if (aIndex + bIndex == 7) {
				x = base.zero();
				y = base.one();
			} else {
				switch (aIndex) {
				case 0: // -6
					switch (bIndex) {
					case 0: // -6
						x = base.getInteger(2);
						y = base.one();
						break;
					case 3: // -1
						x = base.one();
						y = base.one();
						break;
					default:
						throw new ArithmeticException("not split!");
					}
					break;
				case 1: // -3
					switch (bIndex) {
					case 1: // -3
						x = base.one();
						y = base.getInteger(2);
						break;
					case 3: // -1
						x = base.getInteger(2);
						y = base.one();
						break;
					default:
						throw new ArithmeticException("not split!");
					}
					break;
				case 2: // -2
					switch (bIndex) {
					case 6: // 3
						x = base.one();
						y = base.one();
						break;
					case 7: // 6
						x = base.getInteger(2);
						y = base.one();
						break;
					default:
						throw new ArithmeticException("not split!");
					}
					break;
				case 3: // -1
					switch (bIndex) {
					case 0: // -6
						x = base.getInteger(3);
						y = base.one();
						break;
					case 1: // -3
						x = base.getInteger(2);
						y = base.one();
						break;
					case 5: // 2
						x = base.one();
						y = base.one();
						break;
					default:
						throw new ArithmeticException("not split!");
					}
					break;
				case 5: // 2
					switch (bIndex) {
					case 3: // -1
						x = base.one();
						y = base.one();
						break;
					case 5: // 2
						x = base.getInteger(2);
						y = base.one();
						break;
					default:
						throw new ArithmeticException("not split!");
					}
					break;
				case 6: // 3
					switch (bIndex) {
					case 2: // -2
						x = base.one();
						y = base.one();
						break;
					case 7: // 6
						x = base.getInteger(3);
						y = base.one();
						break;
					default:
						throw new ArithmeticException("not split!");
					}
					break;
				case 7: // 6
					switch (bIndex) {
					case 2: // -2
						x = base.getInteger(2);
						y = base.one();
						break;
					case 6: // 3
						x = base.getInteger(3);
						y = base.one();
						break;
					default:
						throw new ArithmeticException("not split!");
					}
					break;
				default:
					throw new ArithmeticException("a is a square!");
				}
			}
		} else {
			if (bIndex == 0) {
				x = base.one();
				y = base.zero();
			} else {
				int s = MiscAlgorithms.jacobiSymbol(BigInteger.valueOf(-1), base.reduction().characteristic());
				if (s == 1) { // -1 is square
					x = base.zero();
					y = base.one();
				} else {
					FiniteField fp2 = reduction.getExtension(reduction.getUnivariatePolynomialRing()
							.getPolynomial(reduction.one(), reduction.zero(), reduction.one())).extension();
					FFE sqrtm1 = fp2.power(fp2.primitiveRoot(),
							fp2.characteristic().subtract(BigInteger.ONE).shiftRight(1));
					x = base.lift(sqrtm1.asPolynomial().univariateCoefficient(0));
					y = base.lift(sqrtm1.asPolynomial().univariateCoefficient(1));
				}
			}
		}
		int maxValueSqrt = Math.max(aValueSqrt, bValueSqrt);
		PAdicField highAccuracyBase = base.withAccuracy(base.getAccuracy() + 2 * (maxValueSqrt));
		if (y.equals(base.zero())) {
			maxValueSqrt = bValueSqrt;
		} else if (x.equals(base.zero())) {
			maxValueSqrt = Math.abs(aValueSqrt - bValueSqrt);
		}
		x = highAccuracyBase.multiply(x, highAccuracyBase.power(base.uniformizer(), maxValueSqrt));
		y = highAccuracyBase.multiply(y, highAccuracyBase.power(base.uniformizer(), maxValueSqrt - aValueSqrt));
		PAdicNumber result = highAccuracyBase.divide(
				highAccuracyBase.subtract(highAccuracyBase.multiply(x, x), highAccuracyBase.multiply(a(), y, y)), b());
		PAdicNumber z = base.round(highAccuracyBase.sqrt(result).keySet().iterator().next(), base.getAccuracy());
		List<Quaternion<PAdicNumber>> candidates = new ArrayList<>();
		candidates.add(normalize(getElement(base.zero(), y, z, base.zero())));
		candidates.add(normalize(getElement(base.zero(), base.negative(y), z, base.zero())));
		candidates.add(normalize(getElement(base.zero(), x, base.zero(), z)));
		candidates.add(normalize(getElement(base.zero(), base.negative(x), base.zero(), z)));
		candidates.add(normalize(getElement(base.zero(), base.zero(), x, y)));
		candidates.add(normalize(getElement(base.zero(), base.zero(), base.negative(x), y)));
		Quaternion<PAdicNumber> splittingElement = null;
		Value minValue = Value.INFINITY;
		for (Quaternion<PAdicNumber> candidate : candidates) {
			Value value = goodSplittingElement(candidate);
			if (splittingElement == null || minValue.compareTo(value) > 0) {
				splittingElement = candidate;
				minValue = value;
			}
		}
		System.err.println(minValue);
		if (minValue.isInfinite()) {
			throw new ArithmeticException("Value too large");
		}
		return splittingElement;
	}
}
