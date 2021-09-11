package fields.interfaces;

import java.math.BigInteger;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeSet;

import fields.integers.Integers.IntE;
import fields.vectors.pivot.PivotStrategy;

public interface Ring<T extends Element<T>> extends MathSet<T> {
	public T zero();

	public T one();

	public BigInteger characteristic();

	public T add(T t1, T t2);

	public T negative(T t);

	public T multiply(T t1, T t2);

	public boolean isUnit(T t);

	public BigInteger getNumberOfUnits();

	public T inverse(T t);

	public boolean isCommutative();

	public boolean isIntegral();

	public boolean isZeroDivisor(T t);

	public boolean isEuclidean();

	public boolean isUniqueFactorizationDomain();

	public static class FactorizationResult<T extends Element<? super T>, U extends Element<? super U>> implements Comparable<FactorizationResult<T, U>> {
		private U unit;
		private SortedMap<T, Integer> factors;

		public FactorizationResult(U unit, SortedMap<T, Integer> factors) {
			this.unit = unit;
			this.factors = factors;
		}

		public U getUnit() {
			return unit;
		}

		public SortedMap<T, Integer> factorMap() {
			return factors;
		}
		
		public Set<T> primeFactors() {
			return factors.keySet();
		}

		public T firstPrimeFactor() {
			return factors.firstKey();
		}

		public T lastPrimeFactor() {
			return factors.lastKey();
		}

		public int multiplicity(T prime) {
			Integer m = factors.get(prime);
			if (m == null) {
				return 0;
			}
			return m.intValue();
		}

		public boolean squareFree() {
			for (T prime : primeFactors()) {
				if (multiplicity(prime) > 1) {
					return false;
				}
			}
			return true;
		}

		public boolean isIrreducible() {
			return factors.size() == 1 && factors.get(factors.firstKey()) == 1;
		}

		public String toString() {
			StringBuilder build = new StringBuilder();
			build.append(unit);
			for (T f : factors.keySet()) {
				build.append("*(" + f + ")");
				int power = factors.get(f);
				if (power != 1) {
					build.append("^" + power);
				}
			}
			return build.toString();
		}

		@Override
		public int compareTo(FactorizationResult<T, U> o) {
			SortedSet<T> combinedPrimes = new TreeSet<>();
			combinedPrimes.addAll(primeFactors());
			combinedPrimes.addAll(o.primeFactors());
			for (T prime : combinedPrimes) {
				int multiplicityHere = multiplicity(prime);
				int multiplicityThere = o.multiplicity(prime);
				if (multiplicityHere != multiplicityThere) {
					return multiplicityHere - multiplicityThere;
				}
			}
			return unit.compareTo(o.unit);
		}

		@SuppressWarnings("unchecked")
		@Override
		public boolean equals(Object obj) {
			if (!(obj instanceof FactorizationResult<?, ?>)) {
				return false;
			}
			return compareTo((FactorizationResult<T, U>) obj) == 0;
		}

	}

	public FactorizationResult<T, T> uniqueFactorization(T t);

	public boolean isIrreducible(T t);

	public boolean isPrime(T t);

	public boolean isSquareFree(T t);

	public List<T> factors(T t);

	public boolean isPrincipalIdealDomain();

	public boolean isDedekindDomain();

	public FactorizationResult<Ideal<T>, Ideal<T>> idealFactorization(Ideal<T> t);

	public boolean isDivisible(T dividend, T divisor);

	public static class QuotientAndRemainderResult<T extends Element<T>> {
		private T quotient;
		private T remainder;

		public QuotientAndRemainderResult(T quotient, T remainder) {
			this.quotient = quotient;
			this.remainder = remainder;
		}

		public T getQuotient() {
			return quotient;
		}

		public T getRemainder() {
			return remainder;
		}

		public String toString() {
			return quotient + " R " + remainder;
		}
	}

	public QuotientAndRemainderResult<T> quotientAndRemainder(T dividend, T divisor);

	public T divideChecked(T dividend, T divisor);

	public T divide(T dividend, T divisor);

	public T remainder(T dividend, T divisor);

	public BigInteger euclidMeasure(T t);

	public PivotStrategy<T> preferredPivotStrategy();

	public Group<T> getAdditiveGroup();

	public Group<T> getMultiplicativeGroup();

	public T add(T t1, T t2, T t3);

	public T add(T t1, T t2, T t3, T t4);

	public T subtract(T minuend, T subtrahend);

	public T getInteger(int n);

	public T getInteger(BigInteger n);

	public T getInteger(IntE n);

	public T multiply(int n, T t);

	public T multiply(BigInteger n, T t);

	public T multiply(IntE n, T t);

	public T multiply(T t1, T t2, T t3);

	public T multiply(int n, T t1, T t2);

	public T multiply(BigInteger n, T t1, T t2);

	public T multiply(IntE n, T t1, T t2);

	public T multiply(int n, T t1, T t2, T t3);

	public T multiply(BigInteger n, T t1, T t2, T t3);

	public T multiply(IntE n, T t1, T t2, T t3);

	public T power(T t, int n);

	public T power(T t, BigInteger n);

	public T projectToUnit(T t);

	public T upToUnit(T t);

	public T numerator(T t);

	public T denominator(T t);

	public T gcd(T t1, T t2);

	public T lcm(T t1, T t2);

	public static class ExtendedEuclideanResult<T extends Element<T>> {
		private T gcd;
		private T coeff1;
		private T coeff2;

		public ExtendedEuclideanResult(T gcd, T coeff1, T coeff2) {
			this.gcd = gcd;
			this.coeff1 = coeff1;
			this.coeff2 = coeff2;
		}

		public T getGcd() {
			return gcd;
		}

		public T getCoeff1() {
			return coeff1;
		}

		public T getCoeff2() {
			return coeff2;
		}

		public String toString() {
			return "GCD: " + gcd + " C1: " + coeff1 + " C2: " + coeff2;
		}
	}

	public ExtendedEuclideanResult<T> extendedEuclidean(T t1, T t2);

	public static class ExtendedEuclideanListResult<T extends Element<T>> {
		private T gcd;
		private List<T> coeffs;

		public ExtendedEuclideanListResult(T gcd, List<T> coeffs) {
			this.gcd = gcd;
			this.coeffs = coeffs;
		}

		public T getGcd() {
			return gcd;
		}

		public List<T> getCoeffs() {
			return coeffs;
		}

		public String toString() {
			return "GCD: " + gcd;
		}
	}

	public ExtendedEuclideanListResult<T> extendedEuclidean(List<T> t);

	public Iterable<T> getUnits();

	public int krullDimension();

	public class IdealResult<T extends Element<T>, I extends Ideal<T>> {
		private List<List<T>> generatorExpressions;
		private List<T> originalGenerators;
		private I ideal;

		public IdealResult(List<List<T>> generatorExpressions, List<T> originalGenerators, I ideal) {
			this.generatorExpressions = generatorExpressions;
			this.originalGenerators = originalGenerators;
			this.ideal = ideal;
		}

		public List<List<T>> getGeneratorExpressions() {
			return generatorExpressions;
		}

		public List<T> getOriginalGenerators() {
			return originalGenerators;
		}

		public I getIdeal() {
			return ideal;
		}

	}
	
	public IdealResult<T, ?> getIdealWithTransforms(List<T> generators);

	@SuppressWarnings("unchecked")
	public IdealResult<T, ?> getIdealWithTransforms(T... generators);

	public Ideal<T> getIdeal(List<T> generators);

	@SuppressWarnings("unchecked")
	public Ideal<T> getIdeal(T... generators);

	public Ideal<T> getUnitIdeal();

	public Ideal<T> getZeroIdeal();

	public Ideal<T> add(Ideal<T> t1, Ideal<T> t2);

	public Ideal<T> multiply(Ideal<T> t1, Ideal<T> t2);

	public Ideal<T> intersect(Ideal<T> t1, Ideal<T> t2);

	public Ideal<T> radical(Ideal<T> t);

	public Ideal<T> power(Ideal<T> t, int power);

	public UnivariatePolynomialRing<T> getUnivariatePolynomialRing();

	public FactorizationResult<Polynomial<T>, T> factorization(UnivariatePolynomial<T> t);

	public boolean hasRoots(Polynomial<T> t);

	public Map<T, Integer> roots(Polynomial<T> t);

	public boolean hasRoot(T t, int n);

	public Map<T, Integer> roots(T t, int n);

	public boolean hasSqrt(T t);

	public Map<T, Integer> sqrt(T t);

	public boolean hasCharacteristicRoot(T t);

	public boolean hasCharacteristicRoot(T t, int power);

	public T characteristicRoot(T t);

	public T characteristicRoot(T t, int power);

	public List<T> adicDevelopment(T t, T base);
}
