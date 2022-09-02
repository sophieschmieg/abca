package fields.interfaces;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeSet;

import fields.integers.Integers.IntE;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import fields.vectors.pivot.PivotStrategy;
import util.Pair;

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

	public boolean isReduced();

	public boolean isIrreducible();

	public boolean isZeroDivisor(T t);

	public boolean isEuclidean();

	public boolean isUniqueFactorizationDomain();

	public class FieldOfFractionsResult<T extends Element<T>, S extends Element<S>> {
		private Ring<T> ring;
		private Field<S> field;
		private MathMap<T, S> embedding;
		private MathMap<S, T> numerator;
		private MathMap<S, T> denominator;
		private MathMap<S, T> asInteger;

		public FieldOfFractionsResult(Ring<T> ring, Field<S> field, MathMap<T, S> embedding, MathMap<S, T> numerator,
				MathMap<S, T> denominator, MathMap<S, T> asInteger) {
			this.ring = ring;
			this.field = field;
			this.embedding = embedding;
			this.numerator = numerator;
			this.denominator = denominator;
			this.asInteger = asInteger;
		}

		public Ring<T> getRing() {
			return ring;
		}

		public Field<S> getField() {
			return field;
		}

		public MathMap<T, S> getEmbedding() {
			return embedding;
		}

		public S embedding(T t) {
			return embedding.evaluate(t);
		}

		public MathMap<S, T> getNumerator() {
			return numerator;
		}

		public T numerator(S t) {
			return numerator.evaluate(t);
		}

		public MathMap<S, T> getDenominator() {
			return denominator;
		}

		public T denominator(S t) {
			return denominator.evaluate(t);
		}

		public MathMap<S, T> getAsInteger() {
			return asInteger;
		}

		public T asInteger(S t) {
			return asInteger.evaluate(t);
		}

		public String toString() {
			return "Q(" + ring + ")";
		}

	}

	public FieldOfFractionsResult<T, ?> fieldOfFractions();

	public class LocalizeResult<T extends Element<T>, S extends Element<S>, U extends Element<U>, R extends Element<R>> {
		private Ring<T> ring;
		private Ideal<T> ideal;
		private LocalRing<S, U, R> localizedRing;
		private MathMap<T, S> embedding;
		private MathMap<S, T> numerator;
		private MathMap<S, T> denominator;
		private MathMap<S, T> asInteger;

		public LocalizeResult(Ring<T> ring, Ideal<T> ideal, LocalRing<S, U, R> localizedRing, MathMap<T, S> embedding,
				MathMap<S, T> numerator, MathMap<S, T> denominator, MathMap<S, T> asInteger) {
			this.ring = ring;
			this.ideal = ideal;
			this.localizedRing = localizedRing;
			this.embedding = embedding;
			this.numerator = numerator;
			this.denominator = denominator;
			this.asInteger = asInteger;
		}

		public Ring<T> getRing() {
			return ring;
		}

		public Ideal<T> getIdeal() {
			return ideal;
		}

		public LocalRing<S, U, R> getLocalizedRing() {
			return localizedRing;
		}

		public MathMap<T, S> getEmbedding() {
			return embedding;
		}

		public S embedding(T t) {
			return embedding.evaluate(t);
		}

		public MathMap<S, T> getNumerator() {
			return numerator;
		}

		public T numerator(S t) {
			return numerator.evaluate(t);
		}

		public MathMap<S, T> getDenominator() {
			return denominator;
		}

		public T denominator(S t) {
			return denominator.evaluate(t);
		}

		public MathMap<S, T> getAsInteger() {
			return asInteger;
		}

		public T asInteger(S t) {
			return asInteger.evaluate(t);
		}

		public String toString() {
			return ideal + "^-1*" + ring;
		}

	}

	public LocalizeResult<T, ?, ?, ?> localizeAtIdeal(Ideal<T> primeIdeal);

	public static class FactorizationResult<T extends Element<? super T>, U extends Element<? super U>>
			implements Comparable<FactorizationResult<T, U>> {
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

		public boolean isPrimePower() {
			return factors.size() == 1;
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

	public static class PrimaryDecompositionResult<T extends Element<T>, I extends Ideal<T>> {
		private List<I> primaryIdeals;
		private List<I> radicals;

		public PrimaryDecompositionResult(List<I> primaryIdeals, List<I> radicals) {
			this.primaryIdeals = primaryIdeals;
			this.radicals = radicals;
		}

		public List<I> getPrimaryIdeals() {
			return primaryIdeals;
		}

		public List<I> getRadicals() {
			return radicals;
		}

		public String toString() {
			return primaryIdeals.toString();
		}
	}

	public PrimaryDecompositionResult<T, ? extends Ideal<T>> primaryDecomposition(Ideal<T> t);

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

	default public T power(T t, IntE n) {
		return power(t, n.getValue());
	}

	public T projectToUnit(T t);

	public T upToUnit(T t);

	DedekindRing<T, ?, ?> asDedekindRing();

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

	public static class BezoutIdentityResult<T extends Element<T>> {
		private List<T> coeff1;
		private List<T> coeff2;

		public BezoutIdentityResult(List<T> coeff1, List<T> coeff2) {
			this.coeff1 = coeff1;
			this.coeff2 = coeff2;
		}

		public List<T> getCoeff1() {
			return coeff1;
		}

		public List<T> getCoeff2() {
			return coeff2;
		}
	}

	/**
	 * For two coprime ideals, computes coefficients c1i and c2j such that sum_i c1i
	 * * g1i + sum_j c2j * g2j = 1
	 * 
	 * @param t1 Ideal of this ring. Coprime to t2.
	 * @param t2 Ideal of this ring. Coprime to t1.
	 * @return two lists of coefficients.
	 */
	public BezoutIdentityResult<T> bezoutIdentity(Ideal<T> t1, Ideal<T> t2);

	public Pair<T, T> bezoutIdentity(T t1, T t2);

	public static class ChineseRemainderPreparation<T extends Element<T>> {
		private List<? extends Ideal<T>> ideals;
		private Ideal<T> product;
		private List<T> multipliers;

		public ChineseRemainderPreparation(List<? extends Ideal<T>> ideals, Ideal<T> product, List<T> multipliers) {
			this.ideals = ideals;
			this.product = product;
			this.multipliers = multipliers;
		}

		public List<? extends Ideal<T>> getIdeals() {
			return ideals;
		}

		public Ideal<T> getProduct() {
			return product;
		}

		public List<T> getMultipliers() {
			return multipliers;
		}

	}

	public ChineseRemainderPreparation<T> prepareChineseRemainderTheorem(List<? extends Ideal<T>> ideals);

	public ChineseRemainderPreparation<T> prepareChineseRemainderTheoremModuli(List<T> moduli);

	public T chineseRemainderTheorem(List<T> elements, ChineseRemainderPreparation<T> preparation);

	public T chineseRemainderTheorem(List<T> elements, List<? extends Ideal<T>> ideals);

	public T chineseRemainderTheoremModuli(List<T> elements, List<T> moduli);

	public boolean coprime(T t1, T t2);

	public boolean coprime(Ideal<T> t1, Ideal<T> t2);

	public Iterable<T> getUnits();

	public int krullDimension();

	public List<? extends Ideal<T>> maximalPrimeIdealChain();

	public List<? extends Ideal<T>> maximalPrimeIdealChain(Ideal<T> start);

	public List<? extends Ideal<T>> maximalPrimeIdealChain(Ideal<T> start, Ideal<T> end);

	public class IdealResult<T extends Element<T>, I extends Ideal<T>> {
		private List<List<T>> generatorExpressions;
		private List<T> originalGenerators;
		private I ideal;
		private List<Vector<T>> syzygies;

		public IdealResult(List<List<T>> generatorExpressions, List<T> originalGenerators, I ideal,
				List<Vector<T>> syzygies) {
			this.generatorExpressions = generatorExpressions;
			this.originalGenerators = originalGenerators;
			this.ideal = ideal;
			this.syzygies = syzygies;
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

		public List<Vector<T>> getSyzygies() {
			return syzygies;
		}

		public String toString() {
			return ideal.toString();
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

	public Ideal<T> getNilRadical();

	public Ideal<T> add(Ideal<T> t1, Ideal<T> t2);

	public Ideal<T> multiply(Ideal<T> t1, Ideal<T> t2);

	public Ideal<T> intersect(Ideal<T> t1, Ideal<T> t2);

	public Ideal<T> radical(Ideal<T> t);

	public Ideal<T> power(Ideal<T> t, int power);

	public <S extends Element<S>> Ideal<T> getIdealEmbedding(Ideal<S> t, MathMap<S, T> map);

	public static class ModuloMaximalIdealResult<T extends Element<T>, S extends Element<S>, R extends Ring<T>, I extends Ideal<T>, F extends Field<S>> {
		private R ring;
		private I ideal;
		private F field;
		private MathMap<T, S> reduction;
		private MathMap<S, T> lift;

		public ModuloMaximalIdealResult(R ring, I ideal, F field, MathMap<T, S> reduction, MathMap<S, T> lift) {
			this.ring = ring;
			this.ideal = ideal;
			this.field = field;
			this.reduction = reduction;
			this.lift = lift;
		}

		public R getRing() {
			return ring;
		}

		public I getIdeal() {
			return ideal;
		}

		public F getField() {
			return field;
		}

		public MathMap<T, S> getReduction() {
			return reduction;
		}

		public S reduce(T t) {
			return reduction.evaluate(t);
		}

		public MathMap<S, T> getLift() {
			return lift;
		}

		public T lift(S t) {
			return lift.evaluate(t);
		}

		@Override
		public String toString() {
			return ring + "/" + ideal;
		}
	}

	public ModuloMaximalIdealResult<T, ?, ?, ?, ?> moduloMaximalIdeal(Ideal<T> ideal);

	public static class ModuloIdealResult<T extends Element<T>, S extends Element<S>> {
		private Ring<T> ring;
		private Ideal<T> ideal;
		private Ring<S> quotientRing;
		private MathMap<T, S> reduction;
		private MathMap<S, T> lift;

		public ModuloIdealResult(Ring<T> ring, Ideal<T> ideal, Ring<S> quotientRing, MathMap<T, S> reduction,
				MathMap<S, T> lift) {
			this.ring = ring;
			this.ideal = ideal;
			this.quotientRing = quotientRing;
			this.reduction = reduction;
			this.lift = lift;
		}

		public Ring<T> getRing() {
			return ring;
		}

		public Ideal<T> getIdeal() {
			return ideal;
		}

		public Ring<S> getQuotientRing() {
			return quotientRing;
		}

		public MathMap<T, S> getReduction() {
			return reduction;
		}

		public S reduce(T t) {
			return reduction.evaluate(t);
		}

		public MathMap<S, T> getLift() {
			return lift;
		}

		public T lift(S t) {
			return lift.evaluate(t);
		}

		@Override
		public String toString() {
			return ring + "/" + ideal;
		}
	}

	public ModuloIdealResult<T, ?> moduloIdeal(Ideal<T> ideal);

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

	public default boolean isIdealMember(List<T> m, T b) {
		return getIdeal(m).contains(b);
	}

	public default Vector<T> asIdealMember(List<T> m, T b) {
		return asSubModuleMember(Matrix.fromRows(Collections.singletonList(new Vector<>(m))), new Vector<>(b));
	}

	public default List<Vector<T>> syzygyProblem(List<T> m) {
		return getIdealWithTransforms(m).getSyzygies();
	}

	public default List<T> simplifyIdealGenerators(List<T> m) {
		List<T> result = new ArrayList<>();
		List<Vector<T>> simplified = simplifySubModuleGenerators(
				Matrix.fromRows(Collections.singletonList(new Vector<>(m))));
		for (Vector<T> generator : simplified) {
			result.add(generator.get(1));
		}
		return result;
	}

	public default boolean isSubModuleMember(Matrix<T> m, Vector<T> b) {
		return isSubModuleMember(m.getModule(this), m, b);
	}

	public default Vector<T> asSubModuleMember(Matrix<T> m, Vector<T> b) {
		return asSubModuleMember(m.getModule(this), m, b);
	}

	public default List<Vector<T>> syzygyProblem(Matrix<T> m) {
		return syzygyProblem(m.getModule(this), m);
	}

	public default List<Vector<T>> simplifySubModuleGenerators(Matrix<T> m) {
		return simplifySubModuleGenerators(m.getModule(this), m);
	}

	public boolean isSubModuleMember(MatrixModule<T> module, Matrix<T> m, Vector<T> b);

	public Vector<T> asSubModuleMember(MatrixModule<T> module, Matrix<T> m, Vector<T> b);

	public List<Vector<T>> syzygyProblem(MatrixModule<T> module, Matrix<T> m);

	public List<Vector<T>> simplifySubModuleGenerators(MatrixModule<T> module, Matrix<T> m);

	public List<T> adicDevelopment(T t, T base);
}
