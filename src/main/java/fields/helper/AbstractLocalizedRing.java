//package fields.helper;
//
//import java.math.BigInteger;
//import java.util.ArrayList;
//import java.util.Collections;
//import java.util.Iterator;
//import java.util.List;
//import java.util.SortedMap;
//import java.util.TreeMap;
//
//import fields.exceptions.InfinityException;
//import fields.interfaces.Element;
//import fields.interfaces.Field;
//import fields.interfaces.Ideal;
//import fields.interfaces.MathMap;
//import fields.interfaces.Polynomial;
//import fields.interfaces.Ring;
//import fields.interfaces.UnivariatePolynomial;
//import fields.interfaces.UnivariatePolynomialRing;
//import util.ConCatMap;
//
//public abstract class AbstractLocalizedRing<T extends Element<T>, S extends Element<S>, R extends Element<R>> extends AbstractRing<S> {
//	private Ring<T> ring;
//	private FieldOfFractionsResult<T, S> fieldOfFractions;
//	private Field<S> field;
//	private Field<R> reduction;
//	private ModuloMaximalIdealResult<T, R> localizedAndModulo;
//	private Ideal<T> ideal;
//	private int dimension;
//
//	public AbstractLocalizedRing(FieldOfFractionsResult<T, S> fieldOfFractions, ModuloMaximalIdealResult<T, R> localizedAndModulo) {
//		this.fieldOfFractions = fieldOfFractions;
//		this.ring = fieldOfFractions.getRing();
//		this.field = fieldOfFractions.getField();
//		this.ideal = localizedAndModulo.getIdeal();
//		this.reduction = localizedAndModulo.getField();
//		this.localizedAndModulo = localizedAndModulo;
//		this.dimension = ring.krullDimension() - ring.maximalPrimeIdealChain(ideal).size() + 1;
//	}
//	
//	@Override
//	public String toString() {
//		return ring.toString() + "_" + ideal.toString();
//	}
//
//	@Override
//	public S zero() {
//		return field.zero();
//	}
//
//	@Override
//	public S one() {
//		return field.one();
//	}
//
//	@Override
//	public BigInteger characteristic() {
//		return ring.characteristic();
//	}
//
//	@Override
//	public S add(S t1, S t2) {
//		return field.add(t1, t2);
//	}
//
//	@Override
//	public S negative(S t) {
//		return field.negative(t);
//	}
//
//	@Override
//	public S multiply(S t1, S t2) {
//		return field.multiply(t1, t2);
//	}
//
//	@Override
//	public boolean isUnit(S t) {
//		return !ideal.contains(fieldOfFractions.getNumerator().evaluate(t));
//	}
//
//	@Override
//	public S inverse(S t) {
//		if (!isUnit(t)) {
//			throw new ArithmeticException("Not invertible");
//		}
//		return field.inverse(t);
//	}
//
//	@Override
//	public boolean isCommutative() {
//		return true;
//	}
//
//	@Override
//	public boolean isIntegral() {
//		return true;
//	}
//
//	@Override
//	public boolean isReduced() {
//		return true;
//	}
//
//	@Override
//	public boolean isIrreducible() {
//		return true;
//	}
//
//	@Override
//	public boolean isZeroDivisor(S t) {
//		return t.equals(zero());
//	}
//
//	@Override
//	public boolean isEuclidean() {
//		return false;
//	}
//
//	public S getEmbedding(T t) {
//		return fieldOfFractions.getEmbedding().evaluate(t);
//	}
//
//	public T getNumerator(S t) {
//		return fieldOfFractions.getNumerator().evaluate(t);
//	}
//
//	public T getDenominator(S t) {
//		return fieldOfFractions.getDenominator().evaluate(t);
//	}
//
//	public T getAsInteger(S t) {
//		return fieldOfFractions.getAsInteger().evaluate(t);
//	}
//
//	public boolean isElement(S t) {
//		return !ideal.contains(fieldOfFractions.getDenominator().evaluate(t));
//	}
//
//	@Override
//	public S gcd(S t1, S t2) {
//		return getEmbedding(ring.gcd(getNumerator(t1), getNumerator(t2)));
//	}
//
//	@Override
//	public ExtendedEuclideanResult<S> extendedEuclidean(S t1, S t2) {
//		ExtendedEuclideanResult<T> ee = ring.extendedEuclidean(getNumerator(t1), getNumerator(t2));
//		return new ExtendedEuclideanResult<>(getEmbedding(ee.getGcd()),
//				getEmbedding(ring.multiply(getDenominator(t1), ee.getCoeff1())),
//				getEmbedding(ring.multiply(getDenominator(t2), ee.getCoeff2())));
//	}
//
//	@Override
//	public boolean isUniqueFactorizationDomain() {
//		return ring.isUniqueFactorizationDomain();
//	}
//
//	@Override
//	public FactorizationResult<S, S> uniqueFactorization(S t) {
//		T numerator = fieldOfFractions.getNumerator().evaluate(t);
//		T denominator = fieldOfFractions.getNumerator().evaluate(t);
//		FactorizationResult<T, T> numeratorFactors = ring.uniqueFactorization(numerator);
//		S unit = inverse(getEmbedding(denominator));
//		SortedMap<S, Integer> factorization = new TreeMap<>();
//		for (T factor : numeratorFactors.primeFactors()) {
//			if (ideal.contains(factor)) {
//				factorization.put(getEmbedding(factor), numeratorFactors.multiplicity(factor));
//			} else {
//				unit = multiply(unit, getEmbedding(factor));
//			}
//		}
//		return new FactorizationResult<>(unit, factorization);
//	}
//
//	@Override
//	public FactorizationResult<Polynomial<S>, S> factorization(UnivariatePolynomial<S> t) {
//		if (!ring.isUniqueFactorizationDomain()) {
//			throw new ArithmeticException("Not a UFD!");
//		}
//		UnivariatePolynomialRing<S> polynomialRing = getUnivariatePolynomialRing();
//		S content = polynomialRing.content(t);
//		t = polynomialRing.contentFree(t);
//		FactorizationResult<S, S> contentFactorization = uniqueFactorization(content);
//		SortedMap<Polynomial<S>, Integer> factors = new TreeMap<>();
//		for (S prime : contentFactorization.primeFactors()) {
//			factors.put(polynomialRing.getEmbedding(prime), contentFactorization.multiplicity(prime));
//		}
//		T lcm = ring.one();
//		for (int i = 0; i <= t.degree(); i++) {
//			lcm = ring.lcm(lcm, fieldOfFractions.getDenominator().evaluate(t.univariateCoefficient(i)));
//		}
//		t = polynomialRing.toUnivariate(polynomialRing.scalarMultiply(getEmbedding(lcm), t));
//		UnivariatePolynomial<T> overRing = ring.getUnivariatePolynomialRing().getEmbedding(t, new MathMap<>() {
//			@Override
//			public T evaluate(S t) {
//				return fieldOfFractions.getAsInteger().evaluate(t);
//			}
//		});
//		S unit = multiply(getEmbedding(ring.getUnivariatePolynomialRing().content(overRing)),
//				contentFactorization.getUnit());
//		FactorizationResult<Polynomial<T>, T> overRingFactors = ring.factorization(overRing);
//		for (Polynomial<T> overRingFactor : overRingFactors.primeFactors()) {
//			int multiplicity = overRingFactors.multiplicity(overRingFactor);
//			Polynomial<S> factor = polynomialRing.getEmbedding(overRingFactor, new MathMap<>() {
//				@Override
//				public S evaluate(T t) {
//					return getEmbedding(t);
//				}
//			});
//			factors.put(factor, multiplicity);
//		}
//		return new FactorizationResult<>(unit, factors);
//	}
//
//	@Override
//	public boolean isPrincipalIdealDomain() {
//		return isDedekindDomain() && isUniqueFactorizationDomain();
//	}
//
//	@Override
//	public boolean isDedekindDomain() {
//		return krullDimension() == 1 /* TODO: && integrallyClosed() */;
//	}
//
//	@Override
//	public boolean isDivisible(S dividend, S divisor) {
//		if (divisor.equals(zero())) {
//			return false;
//		}
//		return isElement(field.divide(dividend, divisor));
//	}
//
//	@Override
//	public QuotientAndRemainderResult<S> quotientAndRemainder(S dividend, S divisor) {
//		if (divisor.equals(zero())) {
//			return new QuotientAndRemainderResult<>(zero(), dividend);
//		}
//		S quotient = field.divide(dividend, divisor);
//		if (isElement(quotient)) {
//			return new QuotientAndRemainderResult<>(quotient, zero());
//		}
//		return new QuotientAndRemainderResult<>(zero(), dividend);
//	}
//
//	@Override
//	public BigInteger euclidMeasure(S t) {
//		return null;
//	}
//
//	@Override
//	public S projectToUnit(S t) {
//		return inverse(getEmbedding(fieldOfFractions.getDenominator().evaluate(t)));
//	}
//
//	@Override
//	public Iterable<S> getUnits() {
//		throw new InfinityException();
//	}
//
//	@Override
//	public int krullDimension() {
//		return dimension;
//	}
//
//	@Override
//	public List<Ideal<S>> maximalPrimeIdealChain(Ideal<S> start) {
//		List<Ideal<T>> chain = ring.maximalPrimeIdealChain(intersectIdealWithRing(start), ideal);
//		List<Ideal<S>> result = new ArrayList<>();
//		for (Ideal<T> chainIdeal : chain) {
//			result.add(getIdeal(chainIdeal));
//		}
//		return result;
//	}
//
//	@Override
//	public List<Ideal<S>> maximalPrimeIdealChain(Ideal<S> start, Ideal<S> end) {
//		Ideal<T> startIdeal = intersectIdealWithRing(start);
//		Ideal<T> endIdeal = intersectIdealWithRing(end);
//		List<Ideal<T>> chain = ring.maximalPrimeIdealChain(startIdeal, endIdeal);
//		List<Ideal<S>> result = new ArrayList<>();
//		for (Ideal<T> chainIdeal : chain) {
//			result.add(getIdeal(chainIdeal));
//		}
//		return result;
//	}
//
//	@Override
//	public LocalizedIdeal<T, S, R> getUnitIdeal() {
//		return getIdeal(ring.getUnitIdeal());
//	}
//
//	@Override
//	public LocalizedIdeal<T, S, R> getZeroIdeal() {
//		return getIdeal(ring.getZeroIdeal());
//	}
//
//	@Override
//	public LocalizedIdeal<T, S, R> getNilRadical() {
//		return getZeroIdeal();
//	}
//
//	public LocalizedIdeal<T, S, R> getIdeal(Ideal<T> ideal) {
//		return new LocalizedIdeal<>(this, ideal);
//	}
//
//	@SuppressWarnings("unchecked")
//	public Ideal<T> intersectIdealWithRing(Ideal<S> ideal) {
//		return ((LocalizedIdeal<T, S, R>) ideal).ideal;
//	}
//
//	@Override
//	public IdealResult<S, LocalizedIdeal<T, S, R>> getIdealWithTransforms(List<S> generators) {
//		List<T> numerators = new ArrayList<>();
//		List<T> denominators = new ArrayList<>();
//		for (S generator : generators) {
//			numerators.add(fieldOfFractions.getNumerator().evaluate(generator));
//			denominators.add(fieldOfFractions.getDenominator().evaluate(generator));
//		}
//		IdealResult<T, ? extends Ideal<T>> numeratorIdeal = ring.getIdealWithTransforms(numerators);
//		if (!ideal.contains(numeratorIdeal.getIdeal())) {
//			for (int i = 0; i < numeratorIdeal.getIdeal().generators().size(); i++) {
//				T generator = numeratorIdeal.getIdeal().generators().get(i);
//				if (!ideal.contains(generator)) {
//					List<T> generatorExpression = numeratorIdeal.getGeneratorExpressions().get(i);
//					List<S> expression = new ArrayList<>();
//					S inverseGenerator = inverse(getEmbedding(generator));
//					for (int j = 0; j < generators.size(); j++) {
//						T e = generatorExpression.get(j);
//						T denom = denominators.get(j);
//						expression.add(multiply(inverseGenerator, getEmbedding(denom), getEmbedding(e)));
//					}
//					return new IdealResult<>(Collections.singletonList(expression), generators, getUnitIdeal());
//				}
//			}
//		}
//		LocalizedIdeal<T, S, R> result = getIdeal(numeratorIdeal.getIdeal());
//		List<List<S>> expressions = new ArrayList<>();
//		for (int i = 0; i < numeratorIdeal.getIdeal().generators().size(); i++) {
//			List<S> expression = new ArrayList<>();
//			for (int j = 0; j < generators.size(); j++) {
//				expression.add(getEmbedding(
//						ring.multiply(denominators.get(j), numeratorIdeal.getGeneratorExpressions().get(i).get(j))));
//			}
//			expressions.add(expression);
//		}
//		return new IdealResult<>(expressions, generators, result);
//	}
//
//	@Override
//	public LocalizedIdeal<T, S, R> intersect(Ideal<S> t1, Ideal<S> t2) {
//		return new LocalizedIdeal<>(this, ring.intersect(intersectIdealWithRing(t1), intersectIdealWithRing(t2)));
//	}
//
//	@Override
//	public Ideal<S> radical(Ideal<S> t) {
//		return new LocalizedIdeal<>(this, ring.radical(intersectIdealWithRing(t)));
//	}
//
//	@Override
//	public Exactness exactness() {
//		return ring.exactness();
//	}
//
//	@Override
//	public S getRandomElement() {
//		S result;
//		do {
//			result = field.getRandomElement();
//		} while (!isElement(result));
//		return result;
//	}
//
//	@Override
//	public boolean isFinite() {
//		return ring.isFinite();
//	}
//
//	@Override
//	public BigInteger getNumberOfElements() throws InfinityException {
//		return ring.getNumberOfElements();
//	}
//
//	@Override
//	public Iterator<S> iterator() {
//		if (ring.isFinite()) {
//			return field.iterator();
//		}
//		throw new InfinityException();
//	}
//	
//	public LocalizedIdeal<T, S, R> maximalIdeal() {
//		return new LocalizedIdeal<>(this, ideal);
//	}
//	
//	@SuppressWarnings("unchecked")
//	public boolean isRegular() {
//		LocalizedIdeal<T, S, R> maximalIdeal = maximalIdeal();
//		LocalizedIdeal<T, S, R> maximalIdealSquared = (LocalizedIdeal<T, S, R>)multiply(maximalIdeal, maximalIdeal);
//	Field<R> red = reduction;
//	
//	}
//
//	
//	@Override
//	public ModuloMaximalIdealResult<S, R> moduloMaximalIdeal(Ideal<S> t) {
//		if (!t.isMaximal()) {
//			throw new ArithmeticException("Not the maximal ideal of the ring!");
//		}
//		return new ModuloMaximalIdealResult<>(this, t, reduction, new MathMap<>() {
//			@Override
//			public R evaluate(S t) {
//				T numerator = fieldOfFractions.getNumerator().evaluate(t);
//				T denominator = fieldOfFractions.getDenominator().evaluate(t);
//				R reducedNumerator = localizedAndModulo.getReduction().evaluate(numerator);
//				R reducedDenominator = localizedAndModulo.getReduction().evaluate(denominator);
//				return reduction.divide(reducedNumerator, reducedDenominator);
//			}
//		}, new ConCatMap<>(localizedAndModulo.getLift(), fieldOfFractions.getEmbedding()));
//	}
//
//	@SuppressWarnings("unchecked")
//	@Override
//	public ModuloMaximalIdealResult<S, ?> localizeAndModOutMaximalIdeal(Ideal<S> ideal) {
//		LocalizedIdeal<T, S, R> localizedIdeal = (LocalizedIdeal<T, S, R>) ideal;
//		return moduloMaximalIdeal(ring.localizeAndModOutMaximalIdeal(localizedIdeal.ideal), localizedIdeal);
//	}
//
//	public static class LocalizedIdeal<T extends Element<T>, S extends Element<S>, R extends Element<R>> extends AbstractIdeal<S> {
//		private AbstractLocalizedRing<T, S, R> localizedRing;
//		private Ideal<T> ideal;
//		private List<S> generators;
//
//		private LocalizedIdeal(AbstractLocalizedRing<T, S, R> localizedRing, Ideal<T> ideal) {
//			super(localizedRing);
//			this.localizedRing = localizedRing;
//			this.ideal = ideal;
//			this.generators = new ArrayList<>();
//			for (T generator : ideal.generators()) {
//				generators.add(localizedRing.getEmbedding(generator));
//			}
//		}
//
//		@Override
//		public boolean isPrimary() {
//			return ideal.isPrimary();
//		}
//
//		@Override
//		public boolean isPrime() {
//			return ideal.isPrime();
//		}
//
//		@Override
//		public boolean isMaximal() {
//			return ideal.equals(localizedRing.ideal);
//		}
//
//		@Override
//		public List<S> generators() {
//			return generators;
//		}
//
//		@Override
//		public List<S> generate(S t) {
//			T numerator = localizedRing.fieldOfFractions.getNumerator().evaluate(t);
//			T denominator = localizedRing.fieldOfFractions.getDenominator().evaluate(t);
//			List<T> generateNumerator = ideal.generate(numerator);
//			List<S> result = new ArrayList<>();
//			S inverseDenominator = localizedRing.inverse(localizedRing.getEmbedding(denominator));
//			for (T generate : generateNumerator) {
//				result.add(localizedRing.multiply(localizedRing.getEmbedding(generate), inverseDenominator));
//			}
//			return result;
//		}
//
//		@Override
//		public S residue(S t) {
//			T numerator = localizedRing.fieldOfFractions.getNumerator().evaluate(t);
//			T denominator = localizedRing.fieldOfFractions.getDenominator().evaluate(t);
//			T residueNumerator = ideal.residue(numerator);
//			return localizedRing.divide(localizedRing.getEmbedding(residueNumerator),
//					localizedRing.getEmbedding(denominator));
//		}
//
//		@Override
//		public boolean contains(S t) {
//			T numerator = localizedRing.fieldOfFractions.getNumerator().evaluate(t);
//			return ideal.contains(numerator);
//		}
//
//		@Override
//		public boolean isFinite() {
//			return localizedRing.isFinite();
//		}
//
//		@Override
//		public BigInteger getNumberOfElements() throws InfinityException {
//			throw new InfinityException();
//		}
//	}
//}
