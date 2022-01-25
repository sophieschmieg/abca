package util;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Random;
import java.util.TreeMap;

import fields.finitefields.PrimeField;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.Ring.QuotientAndRemainderResult;
import fields.local.PAdicField;
import fields.local.PAdicField.PAdicNumber;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.polynomials.Monomial;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;

public class MiscAlgorithms {
	/**
	 * Finds the solutions to a set of modular equations.
	 * 
	 * @param moduloValues the value x the result is supposed to have modulo the
	 *                     corresponding moduli.
	 * @param moduli       the moduli m for the equations.
	 * @param balanced     the solution is unique modulo the lcm of the moduli. If
	 *                     balanced == true, the returned solution will be between
	 *                     -lcm/2 and lcm/2, otherwise between 0 and lcm.
	 * @return n with n = x mod m.
	 */
	public static BigInteger chineseRemainder(List<BigInteger> moduloValues, List<BigInteger> moduli,
			boolean balanced) {
		if (moduloValues.size() != moduli.size()) {
			throw new ArithmeticException("Lists need to be same size for chinese remainder theorem");
		}
		BigInteger modulus = BigInteger.ONE;
		for (BigInteger m : moduli) {
			modulus = modulus.multiply(m).divide(modulus.gcd(m));
		}
		BigInteger result = BigInteger.ZERO;
		for (int i = 0; i < moduli.size(); i++) {
			BigInteger m = moduli.get(i);
			BigInteger multiplier = modulus.divide(m);
			multiplier = multiplier.multiply(multiplier.modInverse(m)).mod(modulus);
			result = result.add(moduloValues.get(i).multiply(multiplier)).mod(modulus);
		}
		if (balanced && result.compareTo(modulus.shiftRight(1)) > 0) {
			result = result.subtract(modulus);
		}
		return result;
	}

	public static int kroneckerSymbol(BigInteger a, BigInteger b) {
		BigInteger zero = BigInteger.ZERO;
		BigInteger one = BigInteger.ONE;
		BigInteger eight = BigInteger.valueOf(8);
		if (b.equals(zero)) {
			if (a.abs().equals(one)) {
				return 1;
			}
			return 0;
		}
		if (b.testBit(0)) {
			int mod8 = a.mod(eight).intValueExact();
			if (mod8 == 1 || mod8 == 7) {
				return kroneckerSymbol(a, b.shiftRight(1));
			} else if (mod8 == 3 || mod8 == 5) {
				return -kroneckerSymbol(a, b.shiftRight(1));
			} else {
				return 0;
			}

		}
		if (b.signum() < 0) {
			if (a.signum() < 0) {
				return -kroneckerSymbol(a, b.negate());
			}
			return kroneckerSymbol(a, b.negate());
		}
		return jacobiSymbol(a, b);
	}

	public static int jacobiSymbol(BigInteger a, BigInteger b) {
		BigInteger zero = BigInteger.ZERO;
		BigInteger one = BigInteger.ONE;
		BigInteger two = BigInteger.TWO;
		BigInteger four = BigInteger.valueOf(4);
		BigInteger eight = BigInteger.valueOf(8);
		if (b.mod(two).equals(zero)) {
			throw new RuntimeException("JacobiSymbol (" + a + "/" + b + ") undefined");
		}
		a = a.mod(b);
		if (a.equals(zero)) {
			return 0;
		}
		// a == -1 mod b
		int aModFour = a.mod(four).intValue();
		int bModFour = b.mod(four).intValue();
		if (a.add(one).equals(b)) {
			if (bModFour == 1) {
				return 1;
			} else if (bModFour == 3) {
				return -1;
			}
		}
		if (a.equals(two)) {
			int bModEight = b.mod(eight).intValue();
			if (bModEight == 1 || bModEight == 7)
				return 1;
			else if (bModEight == 3 || bModEight == 5) {
				return -1;
			}
		}
		if (aModFour % 2 == 0) {
			return jacobiSymbol(two, b) * jacobiSymbol(a.shiftRight(1), b);
		}
		int inv = jacobiSymbol(b, a);
		if (aModFour == 3 && bModFour == 3) {
			return -inv;
		}
		return inv;
	}

	public static String toHexDigit(byte b) {
		String s = Integer.toHexString(Byte.toUnsignedInt(b));
		while (s.length() < 2) {
			s = "0" + s;
		}
		return s;
	}

	public static String toHex(byte[] bytes) {
		StringBuffer buf = new StringBuffer();
		buf.append("0x");
		for (byte b : bytes) {
			buf.append(toHexDigit(b));
		}
		return buf.toString();
	}

	public static BigInteger randomBigInteger(Random rand, BigInteger max) {
		if (max.compareTo(BigInteger.ZERO) <= 0) {
			throw new ArithmeticException("Can't find an integer between 0 and " + max);
		}
		BigInteger num;
		do {
			num = new BigInteger(max.bitLength(), rand);
		} while (num.compareTo(BigInteger.ZERO) < 0 || num.compareTo(max) >= 0);
		return num;
	}

	public static <T extends Element<T>, S extends Element<S>> Polynomial<S> mapPolynomial(Polynomial<T> t,
			MathMap<T, S> map, PolynomialRing<S> codomainRing) {
		Map<Monomial, S> coeff = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			coeff.put(m, map.evaluate(t.coefficient(m)));
		}
		return codomainRing.getPolynomial(coeff);
	}

	public static <T extends Element<T>> MathMap<T, T> identity() {
		return new MathMap<>() {
			@Override
			public T evaluate(T t) {
				return t;
			}
		};
	}

	public static int DivRoundDown(int dividend, int divisor) {
		if (divisor < 0) {
			dividend = -dividend;
			divisor = -divisor;
		} else if (divisor == 0) {
			throw new ArithmeticException("Div by 0");
		}
		if (dividend > 0) {
			return dividend / divisor;
		}
		return (dividend - divisor + 1) / divisor;
	}

	public static int DivRoundUp(int dividend, int divisor) {
		if (divisor < 0) {
			dividend = -dividend;
			divisor = -divisor;
		} else if (divisor == 0) {
			throw new ArithmeticException("Div by 0");
		}
		if (dividend > 0) {
			return (dividend + divisor - 1) / divisor;
		}
		return dividend / divisor;
	}

	public static <T> List<List<T>> crossProduct(List<List<T>> setOfSets) {
		if (setOfSets.isEmpty()) {
			return Collections.singletonList(Collections.emptyList());
		}
		List<T> firstSet = setOfSets.get(0);
		List<List<T>> remainingSets = setOfSets.subList(1, setOfSets.size());
		List<List<T>> product = crossProduct(remainingSets);
		List<List<T>> result = new ArrayList<>();
		for (T element : firstSet) {
			for (List<T> productSet : product) {
				List<T> productElement = new ArrayList<>();
				productElement.add(element);
				productElement.addAll(productSet);
				result.add(productElement);
			}
		}
		return result;
	}

	private static Map<IntE, Map<NFE, Pair<Integer, Integer>>> periods = new TreeMap<>();

	public static Pair<Integer, Integer> continuedFractionPeriod(NumberField nf, NFE sqrt) {
		if (!periods.containsKey(nf.discriminant())) {
			periods.put(nf.discriminant(), new TreeMap<>());
		}
		if (periods.get(nf.discriminant()).containsKey(sqrt)) {
			return periods.get(nf.discriminant()).get(sqrt);
		}
		if (nf.degree() != 2) {
			throw new ArithmeticException("Not a periodic continued fraction!");
		}
		Integers z = Integers.z();
		IntE discriminant = nf.discriminant();
		QuotientAndRemainderResult<IntE> qr = z.quotientAndRemainder(discriminant, z.getInteger(4));
		if (qr.getRemainder().equals(z.zero())) {
			discriminant = qr.getQuotient();
		}
		MathMap<NFE, IntE> roundDown = roundDownSqrt(discriminant);
		NFE state = sqrt;
		Map<NFE, Integer> seenStates = new TreeMap<>();
		int index = 0;
		while (!state.equals(nf.zero())) {
			if (seenStates.containsKey(state)) {
				int head = seenStates.get(state);
				int period = index - head;
				Pair<Integer, Integer> result = new Pair<>(head, period);
				periods.get(nf.discriminant()).put(sqrt, result);
				return result;
			}
			seenStates.put(state, index);
			IntE next = roundDown.evaluate(state);
			state = nf.subtract(state, nf.getInteger(next));
			if (!state.equals(nf.zero())) {
				state = nf.inverse(state);
			}
			index++;
		}
		throw new ArithmeticException("Did not find period?");
	}

	public static <T extends Element<T>> Iterator<IntE> continuedFraction(Field<T> f, T t, MathMap<T, IntE> roundDown) {
		return new Iterator<>() {
			private T state = t;

			@Override
			public boolean hasNext() {
				return !state.equals(f.zero());
			}

			@Override
			public IntE next() {
				IntE next = roundDown.evaluate(state);
				state = f.subtract(state, f.getInteger(next));
				if (!state.equals(f.zero())) {
					state = f.inverse(state);
				}
				return next;
			}
		};
	}

	public static <T extends Element<T>> Iterator<Fraction> continuedFractionApproximation(Field<T> f, T t,
			MathMap<T, IntE> roundDown) {
		Integers z = Integers.z();
		return new Iterator<>() {
			private Iterator<IntE> it = continuedFraction(f, t, roundDown);
			private IntE prevPrevNumerator = z.zero();
			private IntE prevPrevDenominator = z.one();
			private IntE prevNumerator = z.one();
			private IntE prevDenominator = z.zero();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public Fraction next() {
				IntE nextCoeff = it.next();
				IntE nextNumerator = z.add(z.multiply(nextCoeff, prevNumerator), prevPrevNumerator);
				IntE nextDenominator = z.add(z.multiply(nextCoeff, prevDenominator), prevPrevDenominator);
				prevPrevNumerator = prevNumerator;
				prevNumerator = nextNumerator;
				prevPrevDenominator = prevDenominator;
				prevDenominator = nextDenominator;
				return Rationals.q().getFraction(nextNumerator, nextDenominator);
			}
		};
	}

	private static MathMap<NFE, IntE> roundDownSqrt(IntE d) {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		return new MathMap<>() {

			@Override
			public IntE evaluate(NFE t) {
				Fraction rational = t.asPolynomial().univariateCoefficient(0);
				Fraction irrational = t.asPolynomial().univariateCoefficient(1);
				IntE denominator = z.lcm(rational.getDenominator(), irrational.getDenominator());
				IntE rationalNumerator = z.divideChecked(z.multiply(rational.getNumerator(), denominator),
						rational.getDenominator());
				IntE irrationalNumerator = z.divideChecked(z.multiply(irrational.getNumerator(), denominator),
						irrational.getDenominator());
				boolean positive = irrationalNumerator.compareTo(z.zero()) > 0;
				IntE square = z.multiply(irrationalNumerator, irrationalNumerator, d);
				IntE positiveSqrt = z.getInteger(square.getValue().sqrt());
				IntE sqrt = positive ? positiveSqrt : z.subtract(z.negative(positiveSqrt), z.one());
				Fraction result = q.getFraction(z.add(rationalNumerator, sqrt), denominator);
				return result.roundDown();
			}
		};
	}

	public static Iterator<IntE> continuedFractionSqrt(IntE d) {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		if (z.hasSqrt(d)) {
			return q.continuedFraction(q.getInteger(z.sqrt(d).keySet().iterator().next()));
		}
		NumberField nf = new NumberField(
				q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(z.negative(d)), q.zero(), q.one()));
		return continuedFraction(nf, nf.alpha(), roundDownSqrt(d));
	}

	public static Iterator<Fraction> continuedFractionApproximationSqrt(IntE d) {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		if (z.hasSqrt(d)) {
			return q.continuedFractionApproximation(q.getInteger(z.sqrt(d).keySet().iterator().next()));
		}
		NumberField nf = new NumberField(
				q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(z.negative(d)), q.zero(), q.one()));
		return continuedFractionApproximation(nf, nf.alpha(), roundDownSqrt(d));
	}

	/**
	 * Returns (X, Y) such that X^2+Y^2=n, if possible
	 */
	public static Optional<Vector<IntE>> sumOfTwoSquares(IntE n) {
		Integers z = Integers.z();
		if (n.compareTo(z.zero()) < 0) {
			return Optional.empty();
		}
		if (z.isPrime(n)) {
			int mod4 = n.getValue().mod(BigInteger.valueOf(4)).intValueExact();
			switch (mod4) {
			case 0:
				throw new ArithmeticException("Not a prime!");
			case 1:
				PrimeField f = PrimeField.getPrimeField(n);
				IntE root = z
						.getInteger(z.centeredLift(f.sqrt(f.getInteger(-1)).keySet().iterator().next(), n.getValue())
								.getValue().abs());
				IntE cofactor = z.divideChecked(z.add(z.multiply(root, root), z.one()), n);
				if (cofactor.equals(z.one())) {
					return Optional.of(new Vector<>(root, z.one()));
				}
				Vector<IntE> cofactorSquares = sumOfTwoSquares(cofactor).get();
				Vector<IntE> firstRow = new Vector<>(cofactorSquares.get(1), cofactorSquares.get(2));
				Vector<IntE> secondRow = new Vector<>(z.negative(cofactorSquares.get(2)), cofactorSquares.get(1));
				List<Vector<IntE>> rows = new ArrayList<>();
				rows.add(firstRow);
				rows.add(secondRow);
				Vector<IntE> rhs = new Vector<>(root, z.one());
				MatrixAlgebra<IntE> m = new FreeModule<>(z, 2).matrixAlgebra();
				Vector<IntE> result = m.solve(Matrix.fromRows(rows), rhs);
				return Optional.of(result);
			case 2:
				return Optional.of(new Vector<>(z.one(), z.one()));
			case 3:
				return Optional.empty();
			default:
				throw new ArithmeticException("Modulus was wrong");
			}
		}
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(n);
		Vector<IntE> result = new Vector<>(z.one(), z.zero());
		for (IntE prime : factors.primeFactors()) {
			int multiplicity = factors.multiplicity(prime);
			int mod2 = multiplicity % 2;
			Vector<IntE> thisFactor = new Vector<>(z.power(prime, multiplicity / 2), z.zero());
			if (mod2 == 1) {
				Optional<Vector<IntE>> decomposed = sumOfTwoSquares(prime);
				if (decomposed.isEmpty()) {
					return Optional.empty();
				}
				thisFactor = new Vector<>(
						z.add(z.multiply(thisFactor.get(1), decomposed.get().get(1)),
								z.multiply(thisFactor.get(2), decomposed.get().get(2))),
						z.subtract(z.multiply(thisFactor.get(1), decomposed.get().get(2)),
								z.multiply(thisFactor.get(2), decomposed.get().get(1))));
			}
			result = new Vector<>(
					z.add(z.multiply(thisFactor.get(1), result.get(1)), z.multiply(thisFactor.get(2), result.get(2))),
					z.subtract(z.multiply(thisFactor.get(1), result.get(2)),
							z.multiply(thisFactor.get(2), result.get(1))));
		}
		return Optional.of(result);
	}

	/**
	 * Solves X^2 + Y^2 + Z^2 = n
	 */
	public static Optional<Vector<IntE>> sumOfThreeSquares(IntE n) {
		Integers z = Integers.z();
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(n);
		IntE two = z.getInteger(2);
		if (factors.multiplicity(two) % 2 == 0) {
			IntE reduced = z.divideChecked(n, z.power(two, factors.multiplicity(two)));
			IntE remainder8 = z.remainder(reduced, z.getInteger(8));
			if (remainder8.intValueExact() == 7) {
				return Optional.empty();
			}
		}
		if (!factors.squareFree()) {
			IntE sqrt = z.one();
			IntE squareFree = z.one();
			for (IntE factor : factors.primeFactors()) {
				sqrt = z.multiply(sqrt, z.power(factor, MiscAlgorithms.DivRoundDown(factors.multiplicity(factor), 2)));
				squareFree = z.multiply(squareFree, z.power(factor, factors.multiplicity(factor) % 2));
			}
			Vector<IntE> squareFreeSolution = sumOfThreeSquares(squareFree).get();
			return Optional.of(new Vector<>(z.multiply(sqrt, squareFreeSolution.get(1)), z.multiply(sqrt, squareFreeSolution.get(2)), z.multiply(sqrt, squareFreeSolution.get(3))));
		}
		int mod8 = z.remainder(n, z.getInteger(8)).intValueExact();
		return null;
	}

	/**
	 * Solves X^2 - dY^2 = n
	 */
	public static List<Vector<IntE>> pellsEquation(IntE d, IntE n, boolean anySign) {
		Integers z = Integers.z();
		FactorizationResult<IntE, IntE> factorsN = z.uniqueFactorization(n);
		List<List<IntE>> possibleN = new ArrayList<>();
		for (IntE prime : factorsN.primeFactors()) {
			List<IntE> possible = new ArrayList<>();
			int exponent = factorsN.multiplicity(prime);
			int min = exponent % 2;
			for (int i = min; i <= exponent; i += 2) {
				possible.add(z.power(prime, i));
			}
			Collections.reverse(possible);
			possibleN.add(possible);
		}
		List<List<IntE>> crossProduct = MiscAlgorithms.crossProduct(possibleN);
		List<Vector<IntE>> results = new ArrayList<>();
		for (List<IntE> possible : crossProduct) {
			IntE nReduced = factorsN.getUnit();
			for (IntE primePower : possible) {
				nReduced = z.multiply(nReduced, primePower);
			}
			List<Vector<IntE>> result = primitivePellsEquation(d, nReduced, anySign);
			IntE sqrt = z.sqrt(z.divideChecked(n, nReduced)).keySet().iterator().next();
			if (sqrt.compareTo(z.zero()) < 0) {
				sqrt = z.negative(sqrt);
			}
			for (Vector<IntE> r : result) {
				Vector<IntE> correctedResult = new Vector<>(z.multiply(sqrt, r.get(1)), z.multiply(sqrt, r.get(2)));
				results.add(correctedResult);
			}
		}
		return results;
	}

	public static List<Vector<IntE>> primitivePellsEquation(IntE d, IntE n, boolean anySign) {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
//		if (z.hasSqrt(n)) {
//			return Collections.singletonList(new Vector<>(z.sqrt(n).keySet().iterator().next(), z.zero()));
//		}
		if (z.hasSqrt(d)) {
			return null;
		}
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(d);
		for (IntE prime : factors.primeFactors()) {
			PrimeField fp = PrimeField.getPrimeField(prime);
			if (!fp.hasSqrt(fp.getInteger(n)) && !(anySign && fp.hasSqrt(fp.getInteger(z.negative(n))))) {
				return Collections.emptyList();
			}
		}
		FactorizationResult<IntE, IntE> factorsN = z.uniqueFactorization(n);
		List<List<BigInteger>> modSqrts = new ArrayList<>();
		List<BigInteger> mods = new ArrayList<>();
		for (IntE prime : factorsN.primeFactors()) {
			int exponent = factorsN.multiplicity(prime);
			List<BigInteger> sqrts = new ArrayList<>();
			if (prime.equals(z.getInteger(2)) && exponent < 3) {
				if (exponent == 1) {
					sqrts.add(d.getValue().mod(BigInteger.TWO));
				} else if (exponent == 2) {
					int mod4 = d.getValue().mod(BigInteger.valueOf(4)).intValueExact();
					if (mod4 == 0) {
						sqrts.add(BigInteger.ZERO);
					} else if (mod4 == 1) {
						sqrts.add(BigInteger.ONE);
						sqrts.add(BigInteger.valueOf(-1));
					}
				}
			} else {
				PAdicField qp = new PAdicField(prime.getValue(), exponent);
				for (PAdicNumber sqrt : qp.sqrt(qp.round(qp.getInteger(d), exponent)).keySet()) {
					sqrts.add(qp.roundToInteger(sqrt, exponent).getValue());
				}
			}
			if (sqrts.isEmpty()) {
				return Collections.emptyList();
			}
			modSqrts.add(sqrts);
			mods.add(prime.getValue().pow(exponent));
		}
		if (d.compareTo(z.zero()) < 0) {
			return null;
		}
		List<List<BigInteger>> modSqrtsCross = MiscAlgorithms.crossProduct(modSqrts);
		List<Vector<IntE>> results = new ArrayList<>();
		for (List<BigInteger> modSqrt : modSqrtsCross) {
			IntE p0 = z.getInteger(MiscAlgorithms.chineseRemainder(modSqrt, mods, true));
			IntE nAbs = n.compareTo(z.zero()) < 0 ? z.negative(n) : n;
			NumberField nf = new NumberField(
					q.getUnivariatePolynomialRing().getPolynomial(q.negative(q.getInteger(d)), q.zero(), q.one()));
			NFE sqrt = nf.divide(nf.add(nf.getInteger(p0), nf.alpha()), nf.getInteger(nAbs));
			Iterator<Fraction> approximation = continuedFractionApproximation(nf, sqrt, roundDownSqrt(d));
			Pair<Integer, Integer> period = continuedFractionPeriod(nf, sqrt);
			int max = period.getFirst() + 2 * period.getSecond();
			int counter = 0;
			while (approximation.hasNext() && counter <= max) {
				counter++;
				Fraction approx = approximation.next();
				IntE a = approx.getNumerator();
				IntE b = approx.getDenominator();
				IntE g = z.subtract(z.multiply(nAbs, a), z.multiply(p0, b));
				IntE gcd = z.gcd(g, b);
				g = z.divideChecked(g, gcd);
				b = z.divideChecked(b, gcd);
				IntE result = z.subtract(z.multiply(g, g), z.multiply(d, b, b));
				if (result.equals(n) || (anySign && result.equals(z.negative(n)))) {
					results.add(new Vector<>(g, b));
				}
			}
		}
		return results;
	}
}
