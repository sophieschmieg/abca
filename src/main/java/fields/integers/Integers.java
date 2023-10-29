package fields.integers;

import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.ConcurrentModificationException;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Random;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.finitefields.ModuloIntegerRing;
import fields.finitefields.ModuloIntegerRing.ModuloIntegerRingElement;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractIdeal;
import fields.helper.AbstractRing;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.DedekindRing;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.TotalOrder;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.PAdicField;
import fields.local.PAdicField.PAdicNumber;
import fields.local.Value;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.GenericUnivariatePolynomial;
import fields.polynomials.Monomial;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.RealIntegerLattice;
import fields.vectors.Vector;
import util.FunctionMathMap;
import util.Identity;
import util.MiscAlgorithms;
import util.PeekableReader;
import util.SingletonSortedMap;

public class Integers extends AbstractRing<IntE> implements DedekindRing<IntE, Fraction, PFE>, TotalOrder<IntE> {
	private static Integers z = new Integers();
	private Map<IntE, FactorizationResult<IntE, IntE>> factorizatonCache;
	private Map<Polynomial<IntE>, SquareFreeFactorizationResult> polynomialFactorizationCache;
	private IntE zero;
	private IntE one;
	private List<IntE> smallPrimes;
	private Map<Integer, Map<Integer, IntE>> binomialCoefficients;

	private Integers() {
		factorizatonCache = new TreeMap<>();
		polynomialFactorizationCache = new TreeMap<>();
		binomialCoefficients = new TreeMap<>();
		zero = new IntE(0);
		one = new IntE(1);
		smallPrimes = new ArrayList<>();
	}

	public static Integers z() {
		return z;
	}

	public static class IntE extends AbstractElement<IntE> {
		private BigInteger value;

		public IntE(int value) {
			this.value = BigInteger.valueOf(value);
		}

		public IntE(BigInteger value) {
			this.value = value;
		}

		@Override
		public int compareTo(IntE other) {
			return this.value.compareTo(((IntE) other).value);
		}

		public BigInteger getValue() {
			return value;
		}

		public String toString() {
			return this.value.toString();
		}

		public int intValueExact() {
			return value.intValueExact();
		}
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	public Iterator<IntE> primes() {
		return new Iterator<>() {
			private IntE prime = zero();
			private Iterator<IntE> it = smallPrimes.iterator();

			@Override
			public boolean hasNext() {
				return true;
			}

			@Override
			public IntE next() {
				if (it != null && it.hasNext()) {
					try {
						prime = it.next();
					} catch (ConcurrentModificationException e) {
						it = smallPrimes.iterator();
						if (!prime.equals(zero())) {
							while (!it.next().equals(prime))
								;
							prime = it.next();
						}
					}
				} else {
					it = null;
					prime = new IntE(prime.value.nextProbablePrime());
					smallPrimes.add(prime);
				}
				return prime;
			}
		};
	}

	public Iterable<IntE> setOfPrimes() {
		return new Iterable<>() {

			@Override
			public Iterator<IntE> iterator() {
				return primes();
			}
		};
	}

	@Override
	public IntE getRandomElement() {
		return new IntE((int) Math.round(new Random().nextGaussian() * 10.0));
	}

	@Override
	public IntE getRandomElement(IntE max) {
		return new IntE(MiscAlgorithms.randomBigInteger(new Random(), max.value));
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<IntE> iterator() {
		return new Iterator<IntE>() {
			private BigInteger it = BigInteger.ZERO;

			@Override
			public boolean hasNext() {
				return true;
			}

			@Override
			public IntE next() {
				IntE next = new IntE(it);
				if (it.compareTo(BigInteger.ZERO) > 0) {
					it = it.negate();
				} else {
					it = it.negate().add(BigInteger.ONE);
				}
				return next;
			}
		};
	}

	@Override
	public IntE zero() {
		return zero;
	}

	@Override
	public IntE one() {
		return one;
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public IntE add(IntE t1, IntE t2) {
		return new IntE(t1.value.add(t2.value));
	}

	@Override
	public IntE negative(IntE t) {
		return new IntE(t.value.negate());
	}

	@Override
	public IntE multiply(IntE t1, IntE t2) {
		return new IntE(t1.value.multiply(t2.value));
	}

	@Override
	public IntE multiply(IntE t1, IntE t2, IntE t3) {
		return new IntE(t1.value.multiply(t2.value).multiply(t3.value));
	}

	@Override
	public IntE getInteger(BigInteger t) {
		return new IntE(t);
	}

	@Override
	public IntE divide(IntE t1, IntE t2) {
		return new IntE(t1.value.divide(t2.value));
	}

	@Override
	public IntE remainder(IntE t1, IntE t2) {
		if (t2.equals(zero())) {
			return t1;
		}
		return new IntE(t1.value.mod(t2.value));
	}

	@Override
	public boolean isUnit(IntE t) {
		return t.value.equals(BigInteger.ONE) || t.value.equals(BigInteger.ONE.negate());
	}

	@Override
	public BigInteger getNumberOfUnits() {
		return BigInteger.TWO;
	}

	@Override
	public Iterable<IntE> getUnits() {
		List<IntE> list = new ArrayList<>();
		list.add(new IntE(1));
		list.add(new IntE(-1));
		return list;
	}

	@Override
	public IntE inverse(IntE t) {
		return t;
	}

	@Override
	public boolean isCommutative() {
		return true;
	}

	@Override
	public boolean isIntegral() {
		return true;
	}

	@Override
	public boolean isReduced() {
		return true;
	}

	@Override
	public boolean isIrreducible() {
		return true;
	}

	@Override
	public boolean isZeroDivisor(IntE t) {
		return t.value.equals(BigInteger.ZERO);
	}

	@Override
	public boolean isEuclidean() {
		return true;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return true;
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return true;
	}

	@Override
	public boolean isDedekindDomain() {
		return true;
	}

	@Override
	public boolean isDivisible(IntE dividend, IntE divisor) {
		if (dividend.equals(zero())) {
			return true;
		}
		if (divisor.equals(zero())) {
			return false;
		}
		return dividend.value.mod(divisor.value.abs()).equals(BigInteger.ZERO);
	}

	@Override
	public BigInteger euclidMeasure(IntE t) {
		return t.value.abs();
	}

	@Override
	public QuotientAndRemainderResult<IntE> quotientAndRemainder(IntE dividend, IntE divisor) {
		if (divisor.equals(zero())) {
			return new QuotientAndRemainderResult<>(zero(), dividend);
		}
		BigInteger divisorValue = divisor.value.abs();
		BigInteger remainder = dividend.value.mod(divisorValue);
		if (remainder.compareTo(divisorValue.shiftRight(1)) > 0) {
			remainder = remainder.subtract(divisorValue);
		}
		BigInteger quotient = dividend.value.subtract(remainder).divide(divisor.value);
		return new QuotientAndRemainderResult<>(new IntE(quotient), new IntE(remainder));
	}

	@Override
	public IntE gcd(IntE t1, IntE t2) {
		return new IntE(t1.value.gcd(t2.value));
	}

	public boolean likelySmooth(IntE n, IntE max) {
		if (n.equals(zero())) {
			return true;
		}
		PartialFactorizationResult pollardRhoResult = pollardRhoFactorization(n, 1,
				max.getValue().sqrt().intValueExact() + 1);
		return !pollardRhoResult.cofactor.equals(n);
	}

//	public List<Vector<IntE>> polynomialSieve(Polynomial<IntE> f, IntE lowerBoundX, IntE upperBoundX, IntE lowerBoundY,
//			IntE upperBoundY, List<IntE> smoothnessBase) {
//		IntE xDiff = subtract(upperBoundX, lowerBoundX);
//		IntE yDiff = subtract(upperBoundY, lowerBoundY);
//		if (xDiff.compareTo(zero()) <= 0 || yDiff.compareTo(zero()) <= 0) {
//			return Collections.emptyList();
//		}
//		if (xDiff.compareTo(getInteger(Integer.MAX_VALUE)) > 0 || yDiff.compareTo(getInteger(Integer.MAX_VALUE)) > 0) {
//			throw new ArithmeticException("Too many candidates!");
//		}
//		Rationals q = Rationals.q();
//		UnivariatePolynomialRing<IntE> polynomials = getUnivariatePolynomialRing();
//		PolynomialRing<IntE> polynomials2 = AbstractPolynomialRing.getPolynomialRing(this, 2, Monomial.GREVLEX);
//		if (!polynomials2.isHomogeneous(f)) {
//			throw new ArithmeticException("Not homogenous");
//		}
//		UnivariatePolynomial<IntE> dehomogenized = polynomials.getEmbedding(polynomials2.dehomogenize(f, 2),
//				new int[] { 0 });
//		List<IntE> zeroList = new ArrayList<>();
//		zeroList.add(null);
//		zeroList.add(zero());
//		UnivariatePolynomial<IntE> zeroed = polynomials.getEmbedding(polynomials2.partiallyEvaluate(f, zeroList),
//				new int[] { 0 });
//		List<Vector<IntE>> candidates = new ArrayList<>();
//		List<List<IntE>> reconstructed = new ArrayList<>();
//		for (IntE i = lowerBoundX; i.compareTo(upperBoundX) < 0; i = add(i, one())) {
//			List<IntE> reconstructedY = new ArrayList<>();
//			reconstructed.add(reconstructedY);
//			for (IntE j = lowerBoundY; j.compareTo(upperBoundY) < 0; j = add(j, one())) {
//				candidates.add(new Vector<>(i, j));
//				reconstructedY.add(one());
//			}
//		}
//		for (IntE prime : smoothnessBase) {
//			PrimeField reduction = PrimeField.getPrimeField(prime);
//			int power = -1;
//			UnivariatePolynomial<PFE> reduced;
//			IntE primePower = one();
//			do {
//				power++;
//				reduced = reduceUnivariatePolynomial(polynomials.divideScalar(f, primePower), prime);
//				primePower = multiply(prime, primePower);
//			} while (reduced.equals(reduction.getUnivariatePolynomialRing().zero()));
//			primePower = divideChecked(primePower, prime);
//			if (power > 0) {
//				for (int i = 0; i < reconstructed.size(); i++) {
//					reconstructed.set(i, multiply(reconstructed.get(i), primePower));
//				}
//			}
//			Set<IntE> roots = new TreeSet<>();
//			for (PFE root : reduction.roots(reduced).keySet()) {
//				roots.add(lift(reduction.subtract(root, reduce(lowerBound, prime.getValue()))));
//			}
//			while (!roots.isEmpty()) {
//				Set<IntE> nextRoots = new TreeSet<>();
//				primePower = multiply(prime, primePower);
//				for (IntE root : roots) {
//					if (root.compareTo(getInteger(reconstructed.size())) >= 0) {
//						continue;
//					}
//					int increment = primePower.compareTo(z.subtract(upperBound, lowerBound)) < 0
//							? primePower.intValueExact()
//							: reconstructed.size();
//					for (int index = root.intValueExact(); index < reconstructed.size(); index += increment) {
//						reconstructed.set(index, multiply(reconstructed.get(index), prime));
//					}
//					IntE lifted = add(root, lowerBound);
//					IntE result = polynomials.evaluate(f, lifted);
//					result = divideChecked(result, prime);
//					if (result.equals(zero()) && primePower.compareTo(subtract(upperBound, lowerBound)) >= 0) {
//						continue;
//					}
//					UnivariatePolynomial<PFE> nextReduced = reduceUnivariatePolynomial(
//							polynomials
//									.divideScalar(
//											polynomials.substitute(f,
//													Collections.singletonList(
//															polynomials.getPolynomial(lifted, primePower))),
//											primePower),
//							prime);
//					if (nextReduced.equals(reduction.getUnivariatePolynomialRing().zero())) {
//						for (int i = 0; i < prime.intValueExact(); i++) {
//							nextRoots.add(add(root, multiply(i, primePower)));
//						}
//					} else {
//						for (PFE nextRoot : reduction.roots(nextReduced).keySet()) {
//							nextRoots.add(add(root, multiply(lift(nextRoot), primePower)));
//						}
//					}
//				}
//				power++;
//				roots = nextRoots;
//			}
//		}
//		List<Vector<IntE>> result = new ArrayList<>();
//		for (Vector<IntE> candidate : candidates) {
//			IntE eval = polynomials.evaluate(f, candidate);
//			if (eval.compareTo(zero()) < 0) {
//				eval = negative(eval);
//			}
//			if (eval.equals(reconstructed.get(subtract(candidate.get(1), lowerBoundX).intValueExact())
//					.get(subtract(candidate.get(2), lowerBoundY).intValueExact()))) {
//				result.add(candidate);
//			}
//		}
//		return result;
//
//	}

	public List<IntE> polynomialSieveSingleVariable(GenericUnivariatePolynomial<IntE> f, IntE lowerBound,
			IntE upperBound, List<IntE> smoothnessBase) {
		if (lowerBound.compareTo(upperBound) >= 0) {
			return Collections.emptyList();
		}
		if (subtract(upperBound, lowerBound).compareTo(getInteger(Integer.MAX_VALUE)) > 0) {
			throw new ArithmeticException("Too many candidates!");
		}
		UnivariatePolynomialRing<IntE> polynomials = getUnivariatePolynomialRing();
		List<IntE> candidates = new ArrayList<>();
		List<IntE> reconstructed = new ArrayList<>();
		for (IntE i = lowerBound; i.compareTo(upperBound) < 0; i = add(i, one())) {
			candidates.add(i);
			reconstructed.add(one());
		}
		for (IntE prime : smoothnessBase) {
			PrimeField reduction = PrimeField.getPrimeField(prime);
			int power = -1;
			UnivariatePolynomial<PFE> reduced;
			IntE primePower = one();
			do {
				power++;
				reduced = reduceUnivariatePolynomial(polynomials.divideScalar(f, primePower), prime);
				primePower = multiply(prime, primePower);
			} while (reduced.equals(reduction.getUnivariatePolynomialRing().zero()));
			primePower = divideChecked(primePower, prime);
			if (power > 0) {
				for (int i = 0; i < reconstructed.size(); i++) {
					reconstructed.set(i, multiply(reconstructed.get(i), primePower));
				}
			}
			Set<IntE> roots = new TreeSet<>();
			for (PFE root : reduction.roots(reduced).keySet()) {
				roots.add(lift(reduction.subtract(root, reduce(lowerBound, prime.getValue()))));
			}
			while (!roots.isEmpty()) {
				Set<IntE> nextRoots = new TreeSet<>();
				primePower = multiply(prime, primePower);
				for (IntE root : roots) {
					if (root.compareTo(getInteger(reconstructed.size())) >= 0) {
						continue;
					}
					int increment = primePower.compareTo(z.subtract(upperBound, lowerBound)) < 0
							? primePower.intValueExact()
							: reconstructed.size();
					for (int index = root.intValueExact(); index < reconstructed.size(); index += increment) {
						reconstructed.set(index, multiply(reconstructed.get(index), prime));
					}
					IntE lifted = add(root, lowerBound);
					IntE result = polynomials.evaluate(f, lifted);
					result = divideChecked(result, prime);
					if (result.equals(zero()) && primePower.compareTo(subtract(upperBound, lowerBound)) >= 0) {
						continue;
					}
					UnivariatePolynomial<PFE> nextReduced = reduceUnivariatePolynomial(
							polynomials
									.divideScalar(
											polynomials.substitute(f,
													Collections.singletonList(
															polynomials.getPolynomial(lifted, primePower))),
											primePower),
							prime);
					if (nextReduced.equals(reduction.getUnivariatePolynomialRing().zero())) {
						for (int i = 0; i < prime.intValueExact(); i++) {
							nextRoots.add(add(root, multiply(i, primePower)));
						}
					} else {
						for (PFE nextRoot : reduction.roots(nextReduced).keySet()) {
							nextRoots.add(add(root, multiply(lift(nextRoot), primePower)));
						}
					}
				}
				power++;
				roots = nextRoots;
			}
		}
		List<IntE> result = new ArrayList<>();
		for (int i = 0; i < candidates.size(); i++) {
			IntE eval = polynomials.evaluate(f, candidates.get(i));
			if (eval.compareTo(zero()) < 0) {
				eval = negative(eval);
			}
			if (eval.equals(reconstructed.get(i))) {
				result.add(candidates.get(i));
			}
		}
		return result;
	}

	public Optional<FactorizationResult<IntE, IntE>> uniqueFactorizationIfSmooth(IntE n, SortedSet<IntE> primes) {
		if (n.equals(zero())) {
			return Optional.empty();
		}
		if (!likelySmooth(n, primes.last())) {
			return Optional.empty();
		}
		SortedMap<IntE, Integer> factorization = new TreeMap<>();
		for (IntE prime : primes) {
			int multiplicity = 0;
			while (isDivisible(n, prime)) {
				n = divideChecked(n, prime);
				multiplicity++;
			}
			factorization.put(prime, multiplicity);
		}
		if (isUnit(n)) {
			return Optional.of(new FactorizationResult<>(n, factorization));
		}
		return Optional.empty();
	}

	@Override
	public boolean isIrreducible(IntE t) {
		return t.value.isProbablePrime(100);
	}

	@Override
	public FactorizationResult<IntE, IntE> uniqueFactorization(IntE n) {
		if (factorizatonCache.containsKey(n)) {
			return factorizatonCache.get(n);
		}
		IntE t = n;
		IntE unit = one();
		if (n.value.compareTo(BigInteger.ZERO) < 0) {
			unit = negative(unit);
			n = negative(n);
		}
		FactorizationResult<IntE, IntE> result = new FactorizationResult<>(unit, uniqueFactorization(n, false));
		factorizatonCache.put(t, result);
		return result;
	}

	private class PartialFactorizationResult {
		private Map<IntE, Integer> factors;
		private IntE cofactor;

		public PartialFactorizationResult(Map<IntE, Integer> factors, IntE cofactor) {
			this.factors = factors;
			this.cofactor = cofactor;
		}
	}

	private SortedMap<IntE, Integer> uniqueFactorization(IntE n, boolean skipNaive) {
		IntE t = n;
		if (n.value.compareTo(BigInteger.ZERO) < 0) {
			n = negative(n);
		}
		if (n.value.equals(BigInteger.ZERO)) {
			throw new ArithmeticException("0 cannot be factorized!");
		}
		if (n.value.isProbablePrime(100)) {
			return SingletonSortedMap.map(n, 1);
		}
		SortedMap<IntE, Integer> decomposition = new TreeMap<>();
		if (!skipNaive) {
			PartialFactorizationResult naive = naivePrimeDecomposition(n, BigInteger.valueOf(7919));
			n = naive.cofactor;
			decomposition.putAll(naive.factors);
		}
		PartialFactorizationResult pollardRho = pollardRhoFactorization(n, 100, 100000);
		decomposition.putAll(pollardRho.factors);
		n = pollardRho.cofactor;
		if (!n.equals(one())) {
			throw new ArithmeticException("Could not factor " + t + ", left with " + n);
		}
		return decomposition;
	}

	private PartialFactorizationResult pollardRhoFactorization(IntE n, int numberOfTries, int numberOfIterations) {
		if (n.equals(one())) {
			return new PartialFactorizationResult(Collections.emptyMap(), n);
		}
		if (n.value.isProbablePrime(100)) {
			return new PartialFactorizationResult(Collections.singletonMap(n, 1), one());
		}
		BigInteger[] sqrt = n.value.sqrtAndRemainder();
		if (sqrt[1].equals(BigInteger.ZERO)) {
			Map<IntE, Integer> sqrtFactors = uniqueFactorization(getInteger(sqrt[0]), true);
			Map<IntE, Integer> result = new TreeMap<>();
			for (IntE factor : sqrtFactors.keySet()) {
				result.put(factor, 2 * sqrtFactors.get(factor));
			}
			return new PartialFactorizationResult(result, one());
		}
		for (int i = 0; i < numberOfTries; i++) {
			IntE a = getRandomElement(n);
			IntE x = getRandomElement(n);
			IntE y = getRandomElement(n);
			for (int j = 0; j < numberOfIterations; j++) {
				IntE gcd = gcd(n, subtract(x, y));
				if (!gcd.equals(one())) {
					if (gcd.equals(n)) {
						break;
					}
					Map<IntE, Integer> result = new TreeMap<>();
					result.putAll(uniqueFactorization(gcd, true));
					result.putAll(uniqueFactorization(divideChecked(n, gcd), true));
					return new PartialFactorizationResult(result, one());
				}
				x = pollardRhoG(x, n, a);
				y = pollardRhoG(pollardRhoG(y, n, a), n, a);
			}
		}
		return new PartialFactorizationResult(Collections.emptyMap(), n);
	}

	private IntE pollardRhoG(IntE input, IntE n, IntE a) {
		return remainder(add(multiply(input, input), a), n);
	}

	private PartialFactorizationResult naivePrimeDecomposition(IntE n, BigInteger max) {
		if (n.equals(one())) {
			return new PartialFactorizationResult(Collections.emptyMap(), n);
		}
		if (n.value.isProbablePrime(100)) {
			return new PartialFactorizationResult(Collections.singletonMap(n, 1), one());
		}
		Map<IntE, Integer> result = new TreeMap<>();
		IntE zero = zero();
		IntE one = one();

		Iterator<IntE> primes = primes();
		IntE prime = primes.next();
		int power = 0;
		while (prime.value.compareTo(n.value.sqrt().add(one.value)) < 0 && prime.value.compareTo(max) < 0) {
			while (remainder(n, prime).equals(zero)) {
				power++;
				n = divideChecked(n, getInteger(prime));
				result.put(prime, power);
			}
			prime = primes.next();
			power = 0;
		}
		return new PartialFactorizationResult(result, n);
	}

	public IntE eulerToitent(IntE t) {
		FactorizationResult<IntE, IntE> factorization = uniqueFactorization(t);
		IntE result = one();
		for (IntE prime : factorization.primeFactors()) {
			result = multiply(power(prime, factorization.multiplicity(prime) - 1), subtract(prime, one()), result);
		}
		return result;
	}

	@Override
	public IntE projectToUnit(IntE t) {
		int signum = t.value.signum();
		if (signum == 0) {
			return one();
		}
		return getInteger(signum);
	}

	@Override
	public List<Vector<IntE>> simplifySubModuleGenerators(MatrixModule<IntE> module, Matrix<IntE> m) {
		return super.simplifySubModuleGenerators(module, m);
//		return RealIntegerLattice
//				.getIntegerLattice(Reals.r(128), m.rows(), super.simplifySubModuleGenerators(module, m), true)
//				.getModuleGenerators();
//		Rationals q = Rationals.q();
//		RealIntegerLattice result = null;
//		for (Vector<IntE> vector : m.asColumnList()) {
//			if (result == null) {
//				result = new RealIntegerLattice(Reals.r(128), vector.dimension(), Collections.singletonList(vector));
//				continue;
//			}
//			if (result.contains(vector)) {
//				continue;
//			}
//			Vector<Fraction> embeddedVector = Vector.mapVector(q.getEmbeddingMap(), vector);
//			if (!q.isSubModuleMember(result.rationalBaseChange(), embeddedVector)) {
//				List<Vector<IntE>> basis = new ArrayList<>();
//				basis.addAll(result.getBasis());
//				basis.add(vector);
//				result = new RealIntegerLattice(result.getVectorSpace(), result.getFreeModule(), basis);
//				continue;
//			}
//			// Mx = v
//			// (d*I, d*x) = R.D.C
//			// (2 3)
//			// (2 -3) . (5/4 5/6) = (5 0)
//			// (12 0), (0 12), (15 10)
//			Vector<Fraction> asFraction = q.asSubModuleMember(result.rationalBaseChange(), embeddedVector);
//			IntE denominator = one();
//			for (Fraction coeff : asFraction.asList()) {
//				denominator = lcm(coeff.getDenominator(), denominator);
//			}
//			FreeModule<IntE> rankModule = new FreeModule<>(this, result.rank());
//			FiniteRationalVectorSpace rankVectorSpace = new FiniteRationalVectorSpace(result.rank());
//			List<Vector<IntE>> generators = new ArrayList<>();
//			for (Vector<IntE> basisVector : result.getBasis()) {
//				generators.add(rankModule.scalarMultiply(denominator,
//						asSubModuleMember(result.integerBaseChange(), basisVector)));
//			}
//			generators.add(Vector.mapVector(q.getAsIntegerMap(),
//					rankVectorSpace.scalarMultiply(denominator.getValue(), asFraction)));
//			Matrix<IntE> generatorMatrix = Matrix.fromColumns(generators);
//			List<Vector<IntE>> newDiagonalBasis = super.simplifySubModuleGenerators(generatorMatrix.getModule(this), generatorMatrix);
//			List<Vector<IntE>> newBasis = new ArrayList<>();
//			for (Vector<IntE> newBasisVector : newDiagonalBasis) {
//				Vector<IntE> inOldBasis = result.fromVector(newBasisVector);
//				List<IntE> divided = new ArrayList<>();
//				for (IntE coeff : inOldBasis.asList()) {
//					divided.add(divideChecked(coeff, denominator));
//				}
//				newBasis.add(new Vector<>(divided));
//			}
//			result = new RealIntegerLattice(result.getVectorSpace(), result.getFreeModule(), newBasis);
//		}
//		return result.getBasis();
	}

	public static class SmallestIntegerSolutionPreparation {
		private Matrix<IntE> generatorMatrix;
		private MatrixModule<IntE> matrixModule;
		private RealIntegerLattice<Real, Vector<Real>> kernelLattice;
		private FreeModule<IntE> solutionSpace;
		private MathMap<Vector<IntE>, Vector<IntE>> projectionMap;

		private SmallestIntegerSolutionPreparation(List<Vector<IntE>> generators, List<Vector<IntE>> modulus) {
			Integers z = Integers.z();
			int accuracy = 128;
			this.projectionMap = new FunctionMathMap<>(
					(Vector<IntE> t) -> new Vector<>(t.asList().subList(0, generators.size())));
			this.solutionSpace = new FreeModule<>(z, generators.size());
			List<Vector<IntE>> actualGenerators = new ArrayList<>();
			actualGenerators.addAll(generators);
			actualGenerators.addAll(modulus);
			this.generatorMatrix = Matrix.fromColumns(actualGenerators);
			this.matrixModule = generatorMatrix.getModule(Integers.z());
			List<Vector<IntE>> kernelBasis = matrixModule.kernelBasis(generatorMatrix);
			if (kernelBasis.size() == 0) {
				return;
			}
			List<Vector<IntE>> projectedKernelBasis = new ArrayList<>();
			for (Vector<IntE> kernelVector : kernelBasis) {
				projectedKernelBasis.add(projectionMap.evaluate(kernelVector));
			}
			Matrix<IntE> projectedKernelMatrix = Matrix.fromColumns(projectedKernelBasis);
			projectedKernelBasis = projectedKernelMatrix.getModule(Integers.z()).imageBasis(projectedKernelMatrix);
			if (projectedKernelBasis.size() == 0) {
				return;
			}
			Reals r = Reals.r(accuracy);
			this.kernelLattice = RealIntegerLattice.getIntegerLattice(r, generators.size(), projectedKernelBasis);
		}

		private Vector<IntE> smallestKernelVector() {
			return kernelLattice.getModuleGenerators().get(0);
		}

		private List<Vector<IntE>> smallestKernelBasis() {
			return kernelLattice.getModuleGenerators();
		}
	}

	public SmallestIntegerSolutionPreparation prepareSmallestIntegerSolution(List<Vector<IntE>> generators) {
		return prepareSmallestIntegerSolution(generators, Collections.emptyList());
	}

	public SmallestIntegerSolutionPreparation prepareSmallestIntegerSolution(List<Vector<IntE>> generators,
			IntE modulus) {
		int dimension = generators.get(0).dimension();
		FreeModule<IntE> free = new FreeModule<>(this, dimension);
		List<Vector<IntE>> latticeGenerators = new ArrayList<>();
		for (Vector<IntE> basisVector : free.getBasis()) {
			latticeGenerators.add(free.scalarMultiply(modulus, basisVector));
		}
		return prepareSmallestIntegerSolution(generators, latticeGenerators);
	}

	public SmallestIntegerSolutionPreparation prepareSmallestIntegerSolution(List<Vector<IntE>> generators,
			List<Vector<IntE>> modulus) {
		return new SmallestIntegerSolutionPreparation(generators, modulus);
	}

	public SmallestIntegerSolutionPreparation prepareSmallestIntegerSolution(List<Vector<IntE>> generators,
			RealIntegerLattice modulus) {
		return prepareSmallestIntegerSolution(generators, modulus.getModuleGenerators());
	}

	public Vector<IntE> smallestKernelVector(SmallestIntegerSolutionPreparation preparation) {
		return preparation.smallestKernelVector();
	}

	public Vector<IntE> smallestKernelVector(List<Vector<IntE>> generators) {
		return smallestKernelVector(prepareSmallestIntegerSolution(generators));
	}

	public Vector<IntE> smallestKernelVector(List<Vector<IntE>> generators, IntE modulus) {
		return smallestKernelVector(prepareSmallestIntegerSolution(generators, modulus));
	}

	public Vector<IntE> smallestKernelVector(List<Vector<IntE>> generators, List<Vector<IntE>> modulus) {
		return smallestKernelVector(prepareSmallestIntegerSolution(generators, modulus));
	}

	public List<Vector<IntE>> smallestKernelBasis(SmallestIntegerSolutionPreparation preparation) {
		return preparation.smallestKernelBasis();
	}

	public List<Vector<IntE>> smallestKernelBasis(List<Vector<IntE>> generators) {
		return smallestKernelBasis(prepareSmallestIntegerSolution(generators));
	}

	public List<Vector<IntE>> smallestKernelBasis(List<Vector<IntE>> generators, IntE modulus) {
		return smallestKernelBasis(prepareSmallestIntegerSolution(generators, modulus));
	}

	public List<Vector<IntE>> smallestKernelBasis(List<Vector<IntE>> generators, List<Vector<IntE>> modulus) {
		return smallestKernelBasis(prepareSmallestIntegerSolution(generators, modulus));
	}

	public Vector<IntE> smallestIntegerSolution(Vector<IntE> target, SmallestIntegerSolutionPreparation preparation) {
		Vector<IntE> solution = preparation.matrixModule.solve(preparation.generatorMatrix, target);
		if (preparation.kernelLattice == null) {
			return solution;
		}
		Vector<IntE> projectedSolution = preparation.projectionMap.evaluate(solution);
		Vector<IntE> closestKernelVector = preparation.kernelLattice.closestLatticePoint(projectedSolution);
		return preparation.solutionSpace.subtract(solution, closestKernelVector);
	}

	public Vector<IntE> smallestIntegerSolution(List<Vector<IntE>> generators, Vector<IntE> target) {
		return smallestIntegerSolution(target, prepareSmallestIntegerSolution(generators));
	}

	public Vector<IntE> smallestIntegerSolution(List<Vector<IntE>> generators, Vector<IntE> target, IntE modulus) {
		return smallestIntegerSolution(target, prepareSmallestIntegerSolution(generators, modulus));
	}

	public Vector<IntE> smallestIntegerSolution(List<Vector<IntE>> generators, Vector<IntE> target,
			List<Vector<IntE>> modulus) {
		return smallestIntegerSolution(target, prepareSmallestIntegerSolution(generators, modulus));
	}

	@Override
	public int krullDimension() {
		return 1;
	}

	public Fraction round(Fraction t, IntE modulus) {
		IntE denominator = one();
		IntE invertible = t.getDenominator();
		IntE gcd;
		do {
			gcd = gcd(invertible, modulus);
			invertible = divide(invertible, gcd);
			denominator = multiply(denominator, gcd);
		} while (!gcd.equals(one()));
		BigInteger mod = modulus.getValue().multiply(denominator.getValue());
		BigInteger numerator = t.getNumerator().getValue().mod(mod).multiply(invertible.getValue().modInverse(mod))
				.mod(mod);
		return Rationals.q().getFraction(numerator, denominator.getValue());
	}

	public DiscreteValuationRing<Fraction, PFE> localize(BigInteger prime) {
		return Rationals.q().withValuation(prime).ringOfIntegers();
	}

	public DiscreteValuationRing<Fraction, PFE> localize(IntE prime) {
		return localize(prime.value);
	}

	@Override
	public DiscreteValuationRing<Fraction, PFE> localize(Ideal<IntE> prime) {
		return localize(prime.generators().get(0));
	}

	public LocalizedFractions localizeAndQuotient(BigInteger prime) {
		return Rationals.q().withValuation(prime);
	}

	public LocalizedFractions localizeAndQuotient(IntE prime) {
		return localizeAndQuotient(prime.getValue());
	}

	@Override
	public LocalizedFractions localizeAndQuotient(Ideal<IntE> maximalIdeal) {
		return localizeAndQuotient(maximalIdeal.generators().get(0));
	}

	@Override
	public LocalizeResult<IntE, Fraction, Fraction, ?> localizeAtIdeal(Ideal<IntE> primeIdeal) {
		FieldOfFractionsResult<IntE, Fraction> fieldOfFractions = fieldOfFractions();
		if (primeIdeal.equals(getZeroIdeal())) {
			return new LocalizeResult<>(this, primeIdeal, Rationals.q(), fieldOfFractions.getEmbedding(),
					fieldOfFractions.getNumerator(), fieldOfFractions.getDenominator(),
					fieldOfFractions.getAsInteger());
		}
		DiscreteValuationRing<Fraction, PFE> result = localize(primeIdeal);
		return new LocalizeResult<>(this, primeIdeal, result, fieldOfFractions.getEmbedding(),
				fieldOfFractions.getNumerator(), fieldOfFractions.getDenominator(), fieldOfFractions.getAsInteger());
	}

	public Value valuation(IntE t, BigInteger prime) {
		if (t.equals(zero())) {
			return Value.INFINITY;
		}
		int valuation = 0;
		BigInteger n = t.value;
		while (n.mod(prime).equals(BigInteger.ZERO)) {
			valuation++;
			n = n.divide(prime);
		}
		return new Value(valuation);
	}

	public Value valuation(IntE t, Ideal<IntE> maximalIdeal) {
		return valuation(t, maximalIdeal.generators().get(0).value);
	}

	@Override
	public Rationals quotientField() {
		return Rationals.q();
	}

	public IntE archimedeanValue(IntE t) {
		if (t.compareTo(zero) < 0) {
			return negative(t);
		}
		return t;
	}

	@Override
	public boolean hasSqrt(IntE t) {
		if (t.equals(zero())) {
			return true;
		}
		if (t.compareTo(zero()) < 0) {
			return false;
		}
		BigInteger[] sqrt = t.value.sqrtAndRemainder();
		return sqrt[1].equals(BigInteger.ZERO);
	}

	@Override
	public Map<IntE, Integer> sqrt(IntE t) {
		if (t.equals(zero())) {
			return Collections.singletonMap(zero(), 2);
		}
		if (t.compareTo(zero()) < 0) {
			return Collections.emptyMap();
		}
		BigInteger[] sqrt = t.value.sqrtAndRemainder();
		if (!sqrt[1].equals(BigInteger.ZERO)) {
			return Collections.emptyMap();
		}
		Map<IntE, Integer> result = new TreeMap<>();
		result.put(z.getInteger(sqrt[0]), 1);
		result.put(z.negative(z.getInteger(sqrt[0])), 1);
		return result;
	}
	
	public Polynomial<IntE> clearFractionsAndContent(Polynomial<Fraction> polynomial) {
		IntE denominator = one();
		for (Monomial m : polynomial.monomials()) {
			denominator = lcm(denominator, polynomial.coefficient(m).getDenominator());
		}
		PolynomialRing<Fraction> polynomialRing = polynomial.getPolynomialRing();
		PolynomialRing<IntE> intPolynomialRing = AbstractPolynomialRing.getPolynomialRing(this,
				polynomialRing.numberOfVariables(), polynomialRing.getComparator());
		polynomial = polynomialRing.multiply(denominator, polynomial);
		Polynomial<IntE> result = intPolynomialRing.getEmbedding(polynomial, new MathMap<>() {
			@Override
			public IntE evaluate(Fraction t) {
				return t.asInteger();
			}
		});
		return intPolynomialRing.contentFree(result);
	}

	// file:///Users/sophie/Downloads/A_sparse_modular_GCD_algorithm_for_polynomials_ove.pdf
	public CoordinateRingElement<IntE> gcdCoordinateRingElements(CoordinateRing<IntE> coordinateRing,
			CoordinateRingElement<IntE> t1, CoordinateRingElement<IntE> t2) {
		if (coordinateRing.getIdeal().generators().size() > 1) {
			throw new ArithmeticException("Coordinate ring is only allowed to have one equation");
		}
		PolynomialRing<IntE> polynomialRing = coordinateRing.getPolynomialRing();
		int freeVariable = polynomialRing.numberOfVariables() - 1;
		int boundVariable = polynomialRing.numberOfVariables();
		Polynomial<IntE> modulus = coordinateRing.getIdeal().generators().get(0);
		if (modulus.degree(freeVariable) != 0) {
			throw new ArithmeticException("Variable " + freeVariable + " is not free!");
		}
		return null;
	}

	@Override
	public FactorizationResult<Polynomial<IntE>, IntE> factorization(UnivariatePolynomial<IntE> t) {
		UnivariatePolynomialRing<IntE> ring = this.getUnivariatePolynomialRing();
		SortedMap<Polynomial<IntE>, Integer> result = new TreeMap<>();
		IntE content = ring.content(t);
		FactorizationResult<IntE, IntE> contentFactorization = uniqueFactorization(content);
		for (IntE contentFactor : contentFactorization.primeFactors()) {
			result.put(ring.getEmbedding(contentFactor), contentFactorization.multiplicity(contentFactor));
		}
		t = ring.contentFree(t);
		FactorizationResult<Polynomial<IntE>, IntE> squareFreeFactors = ring.squareFreeFactorization(t);
		IntE unit = multiply(contentFactorization.getUnit(), squareFreeFactors.getUnit());
		for (Polynomial<IntE> sff : squareFreeFactors.primeFactors()) {
			int power = squareFreeFactors.multiplicity(sff);
			SquareFreeFactorizationResult factors = factorizeSquareFree(ring.toUnivariate(sff));
			unit = multiply(unit, factors.unit);
			for (Polynomial<IntE> factor : factors.factors) {
				result.put(factor, power);
			}
		}
		return new FactorizationResult<>(unit, result);
	}

	private class CombinedFactors {
		private SortedSet<Integer> usedIndeces;
		private Polynomial<PAdicNumber> combined;

		public CombinedFactors(SortedSet<Integer> usedIndeces, Polynomial<PAdicNumber> combined) {
			this.usedIndeces = usedIndeces;
			this.combined = combined;
		}
	}

	private enum Method {
		LINEAR, QUADRATIC, OKUTSU;
	}

	private class SquareFreeFactorizationResult {
		private List<Polynomial<IntE>> factors;
		private IntE unit;

		private SquareFreeFactorizationResult(List<Polynomial<IntE>> factors, IntE unit) {
			this.factors = factors;
			this.unit = unit;
		}
	}

	private SquareFreeFactorizationResult factorizeSquareFree(UnivariatePolynomial<IntE> t) {
		if (t.degree() == 0) {
			return new SquareFreeFactorizationResult(Collections.singletonList(t), one());
		}
		if (t.degree() < 0) {
			throw new ArithmeticException("Cannot factor 0");
		}
		if (t.degree() == 1) {
			return new SquareFreeFactorizationResult(
					Collections.singletonList(getUnivariatePolynomialRing().contentFree(t)),
					getUnivariatePolynomialRing().content(t));
		}
		if (polynomialFactorizationCache.containsKey(t)) {
			return polynomialFactorizationCache.get(t);
		}
		UnivariatePolynomial<IntE> f = t;
		BigInteger limit = BigInteger.ZERO;
		for (Monomial m : t.monomials()) {
			IntE c = t.coefficient(m);
			BigInteger l = c.getValue().abs();
			if (l.compareTo(limit) > 0) {
				limit = l;
			}
		}
		limit = limit.add(BigInteger.ONE);
		limit = limit.shiftLeft(2);
		for (IntE prime : setOfPrimes()) {
			if (BigInteger.valueOf(t.degree()).mod(prime.value).equals(BigInteger.ZERO)) {
				continue;
			}
			int requiredAccuracy = (int) (Math.log(limit.doubleValue()) / Math.log(prime.value.doubleValue()))
					+ t.degree();
			PAdicField qp = new PAdicField(prime.value, requiredAccuracy);
			DiscreteValuationRing<PAdicNumber, PFE> zp = qp.ringOfIntegers();
			UnivariatePolynomialRing<PAdicNumber> zpr = zp.getUnivariatePolynomialRing();
			UnivariatePolynomial<PAdicNumber> fp = zpr.getEmbedding(t, new MathMap<>() {
				@Override
				public PAdicNumber evaluate(IntE t) {
					return qp.getInteger(t);
				}
			});
			if (!zp.hasGoodReduction(fp)) {
				continue;
			}
			UnivariatePolynomialRing<PFE> reducedRing = zp.reduction().getUnivariatePolynomialRing();
			UnivariatePolynomial<PFE> reduced = zp.reduceUnivariatePolynomial(fp);
			List<Polynomial<PAdicNumber>> liftedFactors = new ArrayList<>();
			// System.err.println("Preparation successful");
			Method method = Method.OKUTSU;
			switch (method) {
			case LINEAR:
				FactorizationResult<Polynomial<PAdicNumber>, PAdicNumber> padicFactorization = zp
						.henselLiftFactorization(zpr.normalize(fp), requiredAccuracy);
				for (Polynomial<PAdicNumber> factor : padicFactorization.primeFactors()) {
					if (factor.degree() > 0) {
						liftedFactors.add(factor);
					}
				}
				break;
			case OKUTSU:
				FactorizationResult<Polynomial<PAdicNumber>, PAdicNumber> padicFactorizationOkutsu = zp
						.factorization(zpr.normalize(fp));
				for (Polynomial<PAdicNumber> factor : padicFactorizationOkutsu.primeFactors()) {
					if (factor.degree() > 0) {
						liftedFactors.add(factor);
					}
				}
				break;
			case QUADRATIC:
				FactorizationResult<Polynomial<PFE>, PFE> factors = zp.reduction().factorization(reduced);
				// System.err.println("Finite Field factorization successful");
				for (Polynomial<PFE> factor : factors.primeFactors()) {
					if (factor.degree() > 0) {
						liftedFactors.add(zp.henselLiftFactor(zpr.normalize(fp), reducedRing.normalize(factor)));
					}
				}
				break;
			}
			// System.err.println("Lifts computed successful");

			Map<Integer, Polynomial<PAdicNumber>> padicFactors = new TreeMap<>();
			for (int i = 0; i < liftedFactors.size(); i++) {
				padicFactors.put(i, liftedFactors.get(i));
			}
			List<CombinedFactors> padicCombinedFactors = new ArrayList<>();
			padicCombinedFactors.add(new CombinedFactors(Collections.emptySortedSet(), zpr.one()));
			List<Polynomial<IntE>> intFactors = new ArrayList<>();
			while (padicFactors.size() > 0) {
				CheckCombinationsResult result = checkCombinations(t, padicFactors, padicCombinedFactors, qp);
				t = result.cofactor;
				intFactors.addAll(result.factors);
				for (int k : result.usedFactors) {
					padicFactors.remove(k);
				}
				padicCombinedFactors = result.combined;
			}
			if (t.degree() != 0) {
				throw new ArithmeticException("Factorization recombination failed!");
			}
			SquareFreeFactorizationResult result = new SquareFreeFactorizationResult(intFactors,
					t.leadingCoefficient());
			polynomialFactorizationCache.put(f, result);
			return result;
//			System.err.println("Fast check not successful");
//
//			liftedFactors.clear();
//			for (int index : padicFactors.keySet()) {
//				liftedFactors.add(padicFactors.get(index));
//				zpr.trace(padicFactors.get(index), 5);
//			}
//			// for (int i = 0; i< liftedFactors.size();i++) {
//			// List<Polynomial<Fraction<IntegerRingElement>>> factors =
//			// }
//
		}
		throw new ArithmeticException("ran out of primes");
	}

	private static class CheckCombinationsResult {
		private List<Polynomial<IntE>> factors = new ArrayList<>();
		private UnivariatePolynomial<IntE> cofactor;
		private Set<Integer> usedFactors = new TreeSet<>();
		private List<CombinedFactors> combined = new ArrayList<>();
	}

	private CheckCombinationsResult checkCombinations(UnivariatePolynomial<IntE> t,
			Map<Integer, Polynomial<PAdicNumber>> padicFactors, List<CombinedFactors> padicCombinedFactors,
			PAdicField qp) {
		CheckCombinationsResult result = new CheckCombinationsResult();
		UnivariatePolynomialRing<PAdicNumber> ring = qp.getUnivariatePolynomialRing();
		result.cofactor = t;
		for (int i : padicFactors.keySet()) {
			Polynomial<PAdicNumber> padicFactor = padicFactors.get(i);
			for (CombinedFactors padicCombinedFactor : padicCombinedFactors) {
				if (padicCombinedFactor.usedIndeces.size() != 0 && padicCombinedFactor.usedIndeces.first() <= i) {
					continue;
				}
				SortedSet<Integer> indeces = new TreeSet<>();
				indeces.addAll(padicCombinedFactor.usedIndeces);
				indeces.add(i);
				Polynomial<PAdicNumber> newCombined = ring.multiply(padicFactor, padicCombinedFactor.combined);
				CheckFactorResult cfr = checkFactor(result.cofactor, newCombined, qp);
				if (cfr.foundFactor) {
					result.factors.add(cfr.factor);
					result.cofactor = cfr.cofactor;
					result.usedFactors.addAll(indeces);
					break;
				} else {
					result.combined.add(new CombinedFactors(indeces, newCombined));
				}
			}
		}
		return result;
	}

	private static class CheckFactorResult {
		private boolean foundFactor = false;
		private Polynomial<IntE> factor = null;
		private UnivariatePolynomial<IntE> cofactor = null;
	}

	private CheckFactorResult checkFactor(UnivariatePolynomial<IntE> t, Polynomial<PAdicNumber> potentialFactor,
			PAdicField qp) {
		CheckFactorResult result = new CheckFactorResult();
		UnivariatePolynomialRing<IntE> ring = getUnivariatePolynomialRing();
		UnivariatePolynomialRing<PAdicNumber> qpRing = qp.getUnivariatePolynomialRing();
		potentialFactor = qpRing.multiply(
				qp.divide(qp.getInteger(t.leadingCoefficient()), potentialFactor.leadingCoefficient()),
				potentialFactor);
		Polynomial<IntE> factor = ring.getEmbedding(potentialFactor, new MathMap<>() {
			@Override
			public IntE evaluate(PAdicNumber number) {
				return qp.roundToInteger(number, qp.getAccuracy());
			}
		});
		factor = ring.contentFree(factor);
		QuotientAndRemainderResult<Polynomial<IntE>> qr = ring.quotientAndRemainder(t, factor);
		if (qr.getRemainder().equals(ring.zero())) {
			result.foundFactor = true;
			result.factor = factor;
			result.cofactor = ring.toUnivariate(qr.getQuotient());
		}
		return result;
	}

	public UnivariatePolynomial<IntE> resolvant(UnivariatePolynomial<IntE> t, Polynomial<IntE> roots) {
		PolynomialRing<IntE> multivariate = AbstractPolynomialRing.getPolynomialRing(this, t.degree(),
				Monomial.GREVLEX);
		return resolvant(t, multivariate.orbit(multivariate.getEmbedding(roots)));
	}

	public UnivariatePolynomial<IntE> resolvant(UnivariatePolynomial<IntE> t, Set<Polynomial<IntE>> orbit) {
		Complex c = Complex.c(128);
		MathMap<IntE, ComplexNumber> embedding = new FunctionMathMap<>((IntE s) -> c.getInteger(s));
		UnivariatePolynomialRing<IntE> polynomialRing = getUnivariatePolynomialRing();
		UnivariatePolynomialRing<ComplexNumber> complexPolynomialRing = c.getUnivariatePolynomialRing();
		PolynomialRing<ComplexNumber> complexMultivariate = AbstractPolynomialRing.getPolynomialRing(c, t.degree(),
				Monomial.GREVLEX);
		Map<ComplexNumber, Integer> complexRoots = c.roots(complexPolynomialRing.getEmbedding(t, embedding));
		List<ComplexNumber> flattened = MiscAlgorithms.flattenMap(complexRoots);
		UnivariatePolynomial<ComplexNumber> complexResult = complexPolynomialRing.one();
		for (Polynomial<IntE> orbitPolynomial : orbit) {
			Polynomial<ComplexNumber> complex = complexMultivariate.getEmbedding(orbitPolynomial, embedding);
			ComplexNumber eval = complexMultivariate.evaluate(complex, flattened);
			complexResult = complexPolynomialRing
					.multiply(complexPolynomialRing.getPolynomial(c.negative(eval), c.one()), complexResult);
		}
		List<IntE> result = new ArrayList<>();
		for (int i = 0; i <= complexResult.degree(); i++) {
			result.add(complexResult.univariateCoefficient(i).realPart().round());
		}
		return polynomialRing.getPolynomial(result);
	}

	public UnivariatePolynomial<IntE> tschirnhausenResolvant(UnivariatePolynomial<IntE> t, Polynomial<IntE> roots) {
		PolynomialRing<IntE> multivariate = AbstractPolynomialRing.getPolynomialRing(this, t.degree(),
				Monomial.GREVLEX);
		return tschirnhausenResolvant(t, multivariate.orbit(multivariate.getEmbedding(roots)));
	}

	public UnivariatePolynomial<IntE> tschirnhausenResolvant(UnivariatePolynomial<IntE> t,
			Set<Polynomial<IntE>> orbit) {

		UnivariatePolynomialRing<IntE> polynomialRing = getUnivariatePolynomialRing();
		UnivariatePolynomial<IntE> resolvant = resolvant(t, orbit);
		if (polynomialRing.isSquareFree(resolvant)) {
			return resolvant;
		}
		for (UnivariatePolynomial<IntE> other : polynomialRing.polynomialSet(t.degree() - 1)) {
			if (other.degree() < 2) {
				continue;
			}
			PolynomialRing<IntE> two = AbstractPolynomialRing.getPolynomialRing(this, 2, Monomial.GREVLEX);
			Polynomial<IntE> t1 = two.getEmbedding(t);
			Polynomial<IntE> t2 = two.subtract(two.getVar(2), two.getEmbedding(other));
			UnivariatePolynomial<IntE> result = polynomialRing.toUnivariate(two.resultant(t1, t2, 1));
			resolvant = resolvant(result, orbit);
			if (polynomialRing.isSquareFree(resolvant)) {
				return resolvant;
			}
		}
		throw new ArithmeticException("Reached end of infinite loop!");
	}

	@Override
	public String toString() {
		return "Z";
	}

	@Override
	public IntE parse(PeekableReader reader) throws IOException {
		StringBuilder build = new StringBuilder();
		boolean first = true;
		while (true) {
			int character = reader.peek();
			if (character < 0) {
				break;
			}
			if (!Character.isDigit((char) character) && !(first && (char) character == '-')) {
				break;
			}
			first = false;
			build.append((char) character);
			reader.skip(1);
		}
		if (build.length() == 0) {
			throw new IOException("No digits found!");
		}
		return getInteger(new BigInteger(build.toString()));
	}

	public IntE binomialCoefficient(int n, int k) {
		if (!binomialCoefficients.containsKey(n)) {
			binomialCoefficients.put(n, new TreeMap<>());
		}
		if (!binomialCoefficients.get(n).containsKey(k)) {
			if (n < 0 || k < 0 || k > n) {
				binomialCoefficients.get(n).put(k, zero());
			} else if (n == 0 && k == 0) {
				binomialCoefficients.get(n).put(k, one());
			} else {
				binomialCoefficients.get(n).put(k,
						add(binomialCoefficient(n - 1, k - 1), binomialCoefficient(n - 1, k)));
			}
		}
		return binomialCoefficients.get(n).get(k);
	}

	@Override
	public IdealResult<IntE, IntegerIdeal> getIdealWithTransforms(List<IntE> generators) {
		if (generators.size() == 0) {
			return new IdealResult<>(Collections.singletonList(Collections.emptyList()), generators,
					new IntegerIdeal(0), Collections.emptyList());
		}
		Matrix<IntE> rowMatrix = Matrix.fromRows(Collections.singletonList(new Vector<>(generators)));
		ExtendedEuclideanListResult<IntE> extendedEuclidean = extendedEuclidean(generators);
		return new IdealResult<>(Collections.singletonList(extendedEuclidean.getCoeffs()), generators,
				new IntegerIdeal(extendedEuclidean.getGcd()), rowMatrix.getModule(this).kernelBasis(rowMatrix));
	}

	@Override
	public IntegerIdeal getUnitIdeal() {
		return getIdeal(Collections.singletonList(one()));
	}

	@Override
	public IntegerIdeal getZeroIdeal() {
		return getIdeal(Collections.emptyList());
	}

	@Override
	public IntegerIdeal getNilRadical() {
		return getZeroIdeal();
	}

	@Override
	public IntegerIdeal getIdeal(List<IntE> generators) {
		return getIdealWithTransforms(generators).getIdeal();
	}

	@Override
	public IntegerIdeal intersect(Ideal<IntE> t1, Ideal<IntE> t2) {
		return getIdeal(Collections.singletonList(lcm(t1.generators().get(0), t2.generators().get(0))));
	}

	@Override
	public IntegerIdeal radical(Ideal<IntE> t) {
		IntE m = t.generators().get(0);
		FactorizationResult<IntE, IntE> factors = uniqueFactorization(m);
		IntE radical = one();
		for (IntE prime : factors.primeFactors()) {
			radical = multiply(radical, prime);
		}
		return getIdeal(Collections.singletonList(radical));
	}

	@Override
	public List<Ideal<IntE>> maximalPrimeIdealChain(Ideal<IntE> start) {
		if (start.equals(getZeroIdeal())) {
			return maximalPrimeIdealChain(start, getIdeal(Collections.singletonList(getInteger(2))));
		}
		return Collections.singletonList(primaryDecomposition(start).getRadicals().get(0));
	}

	@Override
	public List<Ideal<IntE>> maximalPrimeIdealChain(Ideal<IntE> start, Ideal<IntE> end) {
		if (!end.contains(start) || !end.isPrime()) {
			throw new ArithmeticException("Invalid arguments!");
		}
		if (start.equals(getZeroIdeal())) {
			if (end.equals(getZeroIdeal())) {
				return Collections.singletonList(getZeroIdeal());
			}
			List<Ideal<IntE>> result = new ArrayList<>();
			result.add(getZeroIdeal());
			result.add(end);
			return result;
		}
		return Collections.singletonList(end);
	}

	@Override
	public PrimaryDecompositionResult<IntE, IntegerIdeal> primaryDecomposition(Ideal<IntE> t) {
		IntE m = t.generators().get(0);
		FactorizationResult<IntE, IntE> factors = uniqueFactorization(m);
		List<IntegerIdeal> result = new ArrayList<>();
		List<IntegerIdeal> radicals = new ArrayList<>();
		for (IntE prime : factors.primeFactors()) {
			result.add(getIdeal(Collections.singletonList(power(prime, factors.multiplicity(prime)))));
			radicals.add(getIdeal(Collections.singletonList(prime)));
		}
		return new PrimaryDecompositionResult<>(result, radicals);
	}

	@Override
	public ModuloMaximalIdealResult<IntE, PFE, Integers, IntegerIdeal, PrimeField> moduloMaximalIdeal(
			Ideal<IntE> ideal) {
		PrimeField fp = PrimeField.getPrimeField(ideal.generators().get(0).getValue());
		return new ModuloMaximalIdealResult<>(this, (IntegerIdeal) ideal, fp, new MathMap<>() {
			@Override
			public PFE evaluate(IntE t) {
				return reduce(t, ideal);
			}
		}, new MathMap<>() {
			@Override
			public IntE evaluate(PFE t) {
				return lift(t);
			}
		});
	}

	@Override
	public ModuloIdealResult<IntE, ?> moduloIdeal(Ideal<IntE> ideal) {
		if (ideal.equals(getZeroIdeal())) {
			return new ModuloIdealResult<>(this, ideal, this, new Identity<>(), new Identity<>());
		}
		if (ideal.isPrime()) {
			ModuloMaximalIdealResult<IntE, PFE, Integers, IntegerIdeal, PrimeField> mod = moduloMaximalIdeal(ideal);
			return new ModuloIdealResult<>(this, ideal, mod.getField(), mod.getReduction(), mod.getLift());
		}
		ModuloIntegerRing result = new ModuloIntegerRing(ideal.generators().get(0).getValue());
		return new ModuloIdealResult<>(this, ideal, result, new MathMap<>() {

			@Override
			public ModuloIntegerRingElement evaluate(IntE t) {
				return result.reduce(t);
			}
		}, new MathMap<>() {

			@Override
			public IntE evaluate(ModuloIntegerRingElement t) {
				return result.lift(t);
			}
		});
	}

	public IntE lift(ModuloIntegerRingElement t) {
		return new IntE(t.getValue());
	}

	public IntE lift(PFE t) {
		return new IntE(t.getValue());
	}

	public IntE lift(PFE t, Ideal<IntE> maximalIdeal) {
		return lift(t);
	}

	public IntE centeredLift(PFE t, BigInteger prime) {
		IntE lift = lift(t);
		IntE alternative = subtract(lift, new IntE(prime));
		if (lift.value.compareTo(alternative.value.abs()) > 0) {
			return alternative;
		}
		return lift;
	}

	public IntE centeredLift(PFE t, IntE prime) {
		return centeredLift(t, prime.value);
	}

	public IntE centeredLift(PFE t, Ideal<IntE> maximalIdeal) {
		return centeredLift(t, maximalIdeal.generators().get(0).value);
	}

	public Polynomial<IntE> liftPolynomial(Polynomial<PFE> t) {
		return AbstractPolynomialRing.getPolynomialRing(this, t.getPolynomialRing().numberOfVariables(),
				t.getPolynomialRing().getComparator()).getEmbedding(t, new MathMap<>() {
					@Override
					public IntE evaluate(PFE t) {
						return lift(t);
					}
				});
	}

	public UnivariatePolynomial<IntE> liftUnivariatePolynomial(Polynomial<PFE> t) {
		return getUnivariatePolynomialRing().getEmbedding(t, new MathMap<>() {
			@Override
			public IntE evaluate(PFE t) {
				return lift(t);
			}
		});
	}

	public UnivariatePolynomial<IntE> centeredLiftUnivariatePolynomial(Polynomial<PFE> t, BigInteger prime) {
		return getUnivariatePolynomialRing().getEmbedding(t, new MathMap<>() {
			@Override
			public IntE evaluate(PFE t) {
				return centeredLift(t, prime);
			}
		});
	}

	public UnivariatePolynomial<IntE> centeredLiftUnivariatePolynomial(Polynomial<PFE> t, int prime) {
		return centeredLiftUnivariatePolynomial(t, BigInteger.valueOf(prime));
	}

	@Override
	public Field<PFE> reduction(Ideal<IntE> maximalIdeal) {
		return PrimeField.getPrimeField(maximalIdeal.generators().get(0).value);
	}

	@Override
	public FieldOfFractionsResult<IntE, Fraction> fieldOfFractions() {
		Rationals q = Rationals.q();
		return new FieldOfFractionsResult<>(this, q, new MathMap<>() {

			@Override
			public Fraction evaluate(IntE t) {
				return q.getEmbedding(t);
			}
		}, new MathMap<>() {

			@Override
			public IntE evaluate(Fraction t) {
				return t.getNumerator();
			}
		}, new MathMap<>() {

			@Override
			public IntE evaluate(Fraction t) {
				return t.getDenominator();
			}
		}, new MathMap<>() {

			@Override
			public IntE evaluate(Fraction t) {
				return t.asInteger();
			}
		});
	}

	@Override
	public DedekindRing<IntE, Fraction, PFE> asDedekindRing() {
		return this;
	}

	@Override
	public boolean isInteger(Fraction t) {
		return isUnit(t.getDenominator());
	}

	@Override
	public IntE asInteger(Fraction t) {
		return t.asInteger();
	}

	@Override
	public IntE getNumerator(Fraction t) {
		return t.getNumerator();
	}

	@Override
	public IntE getDenominator(Fraction t) {
		return t.getDenominator();
	}

	public PFE reduce(IntE t, BigInteger prime) {
		return PrimeField.getPrimeField(prime).reduce(t);
	}

	public PFE reduce(IntE t, IntE prime) {
		return PrimeField.getPrimeField(prime).reduce(t);
	}

	public PFE reduce(IntE t, Ideal<IntE> prime) {
		return reduce(t, prime.generators().get(0).value);
	}

	public Polynomial<PFE> reducePolynomial(Polynomial<IntE> t, BigInteger prime) {
		return AbstractPolynomialRing.getPolynomialRing(PrimeField.getPrimeField(prime),
				t.getPolynomialRing().numberOfVariables(), t.getPolynomialRing().getComparator())
				.getEmbedding(t, new MathMap<>() {
					@Override
					public PFE evaluate(IntE t) {
						return reduce(t, prime);
					}
				});
	}

	public Polynomial<PFE> reducePolynomial(Polynomial<IntE> t, IntE prime) {
		return reducePolynomial(t, prime.getValue());
	}

	public Polynomial<PFE> reducePolynomial(Polynomial<IntE> t, Ideal<IntE> prime) {
		return reducePolynomial(t, prime.generators().get(0));
	}

	public UnivariatePolynomial<PFE> reduceUnivariatePolynomial(UnivariatePolynomial<IntE> t, BigInteger prime) {
		return PrimeField.getPrimeField(prime).getUnivariatePolynomialRing().getEmbedding(t, new MathMap<>() {
			@Override
			public PFE evaluate(IntE t) {
				return reduce(t, prime);
			}
		});
	}

	public UnivariatePolynomial<PFE> reduceUnivariatePolynomial(UnivariatePolynomial<IntE> t, IntE prime) {
		return reduceUnivariatePolynomial(t, prime.getValue());
	}

	public UnivariatePolynomial<PFE> reduceUnivariatePolynomial(UnivariatePolynomial<IntE> t, Ideal<IntE> prime) {
		return reduceUnivariatePolynomial(t, prime.generators().get(0));
	}

	public boolean isEisenstein(GenericUnivariatePolynomial<IntE> p, BigInteger prime) {
		if (!valuation(p.leadingCoefficient(), prime).equals(new Value(0))) {
			return false;
		}
		Value one = new Value(1);
		if (!valuation(p.univariateCoefficient(0), prime).equals(one)) {
			return false;
		}
		for (int i = 1; i < p.degree(); i++) {
			if (valuation(p.univariateCoefficient(i), prime).compareTo(one) < 0) {
				return false;
			}
		}
		return true;
	}

	public static class IntegerIdeal extends AbstractIdeal<IntE> {
		private IntE m;

		public IntegerIdeal(IntE m) {
			super(z);
			this.m = z.archimedeanValue(m);
		}

		public IntegerIdeal(int m) {
			this(BigInteger.valueOf(m));
		}

		public IntegerIdeal(BigInteger m) {
			this(new IntE(m));
		}

		@Override
		public boolean isFinite() {
			return false;
		}

		@Override
		public BigInteger getNumberOfElements() throws InfinityException {
			throw new InfinityException();
		}

//		@Override
//		public List<List<IntE>> nonTrivialCombinations(List<IntE> s) {
//			Matrix<IntE> m = new Matrix<>(Collections.singletonList(s));
//			MatrixModule<IntE> mm = new MatrixModule<>(z(), 1, s.size());
//			List<Vector<IntE>> kernelBasis = mm.kernelBasis(m);
//			List<List<IntE>> result = new ArrayList<>();
//			for (Vector<IntE> basisVector : kernelBasis) {
//				result.add(basisVector.asList());
//			}
//			return result;
//		}

		@Override
		public List<IntE> generators() {
			return Collections.singletonList(m);
		}

		@Override
		public List<IntE> generate(IntE t) {
			if (m.equals(zero())) {
				return Collections.singletonList(zero());
			}
			return Collections.singletonList(z.getInteger(t.value.divide(m.value)));
		}

		@Override
		public IntE residue(IntE t) {
			if (m.equals(zero())) {
				return t;
			}
			return new IntE(t.value.mod(m.value));
		}

		@Override
		public boolean isPrime() {
			return m.equals(z.zero()) || isMaximal();
		}

		@Override
		public boolean isMaximal() {
			return m.value.isProbablePrime(100);
		}

		@Override
		public boolean isPrimary() {
			return Integers.z().uniqueFactorization(m).primeFactors().size() == 1;
		}

		@Override
		public boolean contains(IntE t) {
			if (m.equals(zero())) {
				return t.equals(zero());
			}
			return t.value.mod(m.value).equals(BigInteger.ZERO);
		}

		@Override
		public Value maximumPowerContains(IntE t) {
			if (t.equals(zero()) || m.equals(z.one())) {
				return Value.INFINITY;
			}
			if (m.equals(zero())) {
				return Value.ZERO;
			}
			FactorizationResult<Ideal<IntE>, Ideal<IntE>> factors = z.idealFactorization(this);
			Value value = Value.INFINITY;
			for (Ideal<IntE> factor : factors.primeFactors()) {
				value = value.min(new Value(z.valuation(t, factor).value() / z.valuation(m, factor).value()));
			}
			return value;
		}
	}
}
