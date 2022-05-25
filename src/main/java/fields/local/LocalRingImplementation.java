package fields.local;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractIdeal;
import fields.helper.AbstractModule;
import fields.helper.AbstractRing;
import fields.helper.FieldEmbedding;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DedekindRing;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationField.OtherVersion;
import fields.interfaces.DiscreteValuationField.Valuation;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Field.Extension;
import fields.interfaces.FieldExtension;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Module;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.interfaces.ValueField.AbsoluteValue;
import fields.polynomials.AbstractPolynomialRing;
import fields.vectors.FreeModule;
import fields.vectors.Vector;
import fields.vectors.pivot.PivotStrategy;
import fields.vectors.pivot.ValuationPivotStrategy;
import util.ConCatMap;
import util.ConstantMap;
import util.FlattenProductVectorMap;
import util.Identity;
import util.MiscAlgorithms;
import util.Pair;
import util.SingletonSortedMap;

public class LocalRingImplementation<T extends Element<T>, S extends Element<S>> extends AbstractRing<T>
		implements DiscreteValuationRing<T, S> {
	private DiscreteValuationField<T, S> field;
	private String asString;

	public LocalRingImplementation(DiscreteValuationField<T, S> field, String asString) {
		this.field = field;
		this.asString = asString;
	}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	@Override
	public boolean isCommutative() {
		return true;
	}

	@Override
	public FieldOfFractionsResult<T, T> fieldOfFractions() {
		return new FieldOfFractionsResult<>(this, localField(), new Identity<>(), new MathMap<>() {

			@Override
			public T evaluate(T t) {
				Value v = valuation(t);
				if (v.compareTo(Value.ZERO) >= 0) {
					return t;
				}
				return localField().multiply(t, power(uniformizer(), -v.value()));
			}
		}, new MathMap<>() {

			@Override
			public T evaluate(T t) {
				Value v = valuation(t);
				if (v.compareTo(Value.ZERO) >= 0) {
					return one();
				}
				return power(uniformizer(), -v.value());
			}
		}, new MathMap<>() {

			@Override
			public T evaluate(T t) {
				Value v = valuation(t);
				if (v.compareTo(Value.ZERO) >= 0) {
					return t;
				}
				throw new ArithmeticException("Not an integer");
			}
		});
	}

	@SuppressWarnings("unchecked")
	@Override
	public LocalizeResult<T, T, T, ?> localizeAtIdeal(Ideal<T> primeIdeal) {
		LocalIdeal<T, S> ideal = (LocalIdeal<T, S>) primeIdeal;
		if (ideal.valuation.equals(Value.ZERO)) {
			FieldOfFractionsResult<T, T> fieldOfFractions = fieldOfFractions();
			return new LocalizeResult<T, T, T, T>(this, ideal, field, fieldOfFractions.getEmbedding(),
					fieldOfFractions.getNumerator(), fieldOfFractions.getDenominator(),
					fieldOfFractions.getAsInteger());
		}
		return new LocalizeResult<>(this, ideal, this, new Identity<>(), new Identity<>(), new ConstantMap<>(one()),
				new Identity<>());
	}

	@Override
	public boolean isComplete() {
		return field.isComplete();
	}

	@Override
	public Real value(T t) {
		return field.value(t);
	}

	@Override
	public AbsoluteValue<T> value() {
		return field.value();
	}

	@Override
	public Value valuation(T t) {
		return field.valuation(t);
	}

	@Override
	public Valuation<T> valuation() {
		return field.valuation();
	}

	@Override
	public Value valuationOfUnivariatePolynomial(Polynomial<T> t) {
		UnivariatePolynomial<T> p = field.getUnivariatePolynomialRing().toUnivariate(t);
		Value value = Value.INFINITY;
		for (int i = 0; i <= t.degree(); i++) {
			value = value.min(valuation(p.univariateCoefficient(i)));
		}
		return value;
	}

	@Override
	public T zero() {
		return field.zero();
	}

	@Override
	public T one() {
		return field.one();
	}

	@Override
	public BigInteger characteristic() {
		return field.characteristic();
	}

	@Override
	public T add(T t1, T t2) {
		return field.add(t1, t2);
	}

	@Override
	public T negative(T t) {
		return field.negative(t);
	}

	@Override
	public T multiply(T t1, T t2) {
		return field.multiply(t1, t2);
	}

	@Override
	public boolean isUnit(T t) {
		return field.valuation(t).equals(new Value(0));
	}

	@Override
	public T inverse(T t) {
		if (!isUnit(t)) {
			throw new ArithmeticException("not invertible in the ring");
		}
		return field.inverse(t);
	}

	@Override
	public T divide(T dividend, T divisor) {
		if (!isDivisible(dividend, divisor)) {
			throw new ArithmeticException("not divisible in the ring");
		}
		return field.divide(dividend, divisor);
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
	public boolean isZeroDivisor(T t) {
		return false;
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
	public FactorizationResult<T, T> uniqueFactorization(T t) {
		int valuation = valuation(t).value();
		if (valuation == 0) {
			return new FactorizationResult<>(t, Collections.emptySortedMap());
		}
		t = divideChecked(t, power(uniformizer(), valuation));
		return new FactorizationResult<>(t, SingletonSortedMap.map(uniformizer(), valuation));
	}

	@Override
	public boolean isDivisible(T dividend, T divisor) {
		if (divisor.equals(zero())) {
			return false;
		}
		return valuation(dividend).compareTo(valuation(divisor)) >= 0;
	}

	@Override
	public PivotStrategy<T> preferredPivotStrategy() {
		return new ValuationPivotStrategy<>(this.valuation());
	}

	@Override
	public QuotientAndRemainderResult<T> quotientAndRemainder(T dividend, T divisor) {
		if (!isDivisible(dividend, divisor)) {
			return new QuotientAndRemainderResult<>(zero(), dividend);
		}
		return new QuotientAndRemainderResult<>(divide(dividend, divisor), zero());
	}

	@Override
	public T gcd(T t1, T t2) {
		if (t1.equals(zero())) {
			return t2;
		}
		if (t2.equals(zero())) {
			return t1;
		}
		return power(uniformizer(), Math.min(valuation(t1).value(), valuation(t2).value()));
	}

	@Override
	public ExtendedEuclideanResult<T> extendedEuclidean(T t1, T t2) {
		T gcd = gcd(t1, t2);
		return new ExtendedEuclideanResult<>(gcd,
				valuation(t1).compareTo(valuation(t2)) <= 0 ? field.divide(gcd, t1) : zero(),
				valuation(t1).compareTo(valuation(t2)) <= 0 ? zero() : field.divide(gcd, t2));
	}

	@Override
	public BigInteger euclidMeasure(T t) {
		if (t.equals(zero())) {
			return null;
		}
		return BigInteger.valueOf(valuation(t).value());
	}

	@Override
	public T projectToUnit(T t) {
		if (t.equals(zero())) {
			return field.one();
		}
		return field.divide(t, power(uniformizer(), valuation(t).value()));
	}

	@Override
	public Iterable<T> getUnits() {
		throw new InfinityException();
	}

	@Override
	public int krullDimension() {
		return 1;
	}

	@Override
	public List<Ideal<T>> maximalPrimeIdealChain(Ideal<T> start) {
		return maximalPrimeIdealChain(start, getIdeal(Collections.singletonList(uniformizer())));
	}

	@Override
	public List<Ideal<T>> maximalPrimeIdealChain(Ideal<T> start, Ideal<T> end) {
		if (!end.contains(start) || !end.isPrime()) {
			throw new ArithmeticException("Invalid arguments");
		}
		if (start.equals(getZeroIdeal())) {
			if (end.equals(getZeroIdeal())) {
				return Collections.singletonList(getZeroIdeal());
			}
			List<Ideal<T>> result = new ArrayList<>();
			result.add(start);
			result.add(end);
			return result;
		}
		return Collections.singletonList(end);
	}

	@Override
	public IdealResult<T, LocalIdeal<T, S>> getIdealWithTransforms(List<T> generators) {
		if (generators.size() == 0) {
			return new IdealResult<>(Collections.singletonList(Collections.emptyList()), generators,
					new LocalIdeal<>(this, Value.INFINITY), Collections.emptyList());
		}
		Value minValue = Value.INFINITY;
		int minValueIndex = -1;
		for (int i = 0; i < generators.size(); i++) {
			Value value = valuation(generators.get(i));
			if (value.compareTo(minValue) < 0) {
				minValue = value;
				minValueIndex = i;
			}
		}
		FreeModule<T> freeModule = new FreeModule<>(this, generators.size());
		if (minValue.isInfinite()) {
			List<T> zeroes = new ArrayList<>();
			for (int i = 0; i < generators.size(); i++) {
				zeroes.add(zero());
			}
			return new IdealResult<>(Collections.singletonList(zeroes), generators,
					new LocalIdeal<>(this, Value.INFINITY), freeModule.getBasis());
		}
		T gcd = power(uniformizer(), minValue.value());
		List<Vector<T>> syzygies = new ArrayList<>();
		List<T> expression = new ArrayList<>();
		T generatorCoefficient = divideChecked(gcd, generators.get(minValueIndex));
		for (int i = 0; i < generators.size(); i++) {
			if (i == minValueIndex) {
				expression.add(generatorCoefficient);
			} else {
				expression.add(zero());
				syzygies.add(freeModule.subtract(freeModule.getUnitVector(i + 1),
						freeModule.scalarMultiply(divideChecked(generators.get(i), generators.get(minValueIndex)),
								freeModule.getUnitVector(minValueIndex + 1))));
			}
		}
		return new IdealResult<>(Collections.singletonList(expression), generators, new LocalIdeal<>(this, minValue),
				syzygies);
	}

	@Override
	public LocalIdeal<T, S> getIdeal(List<T> generators) {
		return getIdealWithTransforms(generators).getIdeal();
	}

	@Override
	public LocalIdeal<T, S> getZeroIdeal() {
		return getIdeal(Collections.emptyList());
	}

	@Override
	public LocalIdeal<T, S> getNilRadical() {
		return getZeroIdeal();
	}

	@Override
	public LocalIdeal<T, S> maximalIdeal() {
		return getIdeal(Collections.singletonList(uniformizer()));
	}

	@SuppressWarnings("unchecked")
	@Override
	public LocalIdeal<T, S> intersect(Ideal<T> t1, Ideal<T> t2) {
		LocalIdeal<T, S> ideal1 = (LocalIdeal<T, S>) t1;
		LocalIdeal<T, S> ideal2 = (LocalIdeal<T, S>) t2;
		return ideal1.valuation.compareTo(ideal2.valuation) >= 0 ? ideal1 : ideal2;
	}

	@Override
	public LocalIdeal<T, S> radical(Ideal<T> t) {
		return getIdeal(Collections.singletonList(uniformizer()));
	}

	@SuppressWarnings("unchecked")
	@Override
	public PrimaryDecompositionResult<T, LocalIdeal<T, S>> primaryDecomposition(Ideal<T> t) {
		return new PrimaryDecompositionResult<>(Collections.singletonList((LocalIdeal<T, S>) t),
				Collections.singletonList(radical(t)));
	}

	@Override
	public ModuloMaximalIdealResult<T, S, LocalRingImplementation<T, S>, Ideal<T>, Field<S>> moduloMaximalIdeal(
			Ideal<T> ideal) {
		if (!ideal.isMaximal()) {
			throw new ArithmeticException("not a maximal ideal!");
		}
		return new ModuloMaximalIdealResult<>(this, ideal, reduction(), new MathMap<>() {
			@Override
			public S evaluate(T t) {
				return reduce(t);
			}
		}, new MathMap<>() {
			@Override
			public T evaluate(S t) {
				return lift(t);
			}
		});
	}

	@Override
	public ModuloIdealResult<T, ?> moduloIdeal(Ideal<T> ideal) {
		throw new UnsupportedOperationException("Not yet implemented!");
	}

	@Override
	public T getRandomElement() {
		return field.getRandomInteger();
	}

	@Override
	public DedekindRing<T, T, S> asDedekindRing() {
		return this;
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
	public Iterator<T> iterator() {
		throw new InfinityException();
	}

	@Override
	public T uniformizer() {
		return field.uniformizer();
	}

	@Override
	public DiscreteValuationField<T, S> localField() {
		return field;
	}

	public boolean isElement(T t) {
		return field.isInteger(t);
	}

	@Override
	public Field<S> reduction() {
		return field.residueField();
	}

	@Override
	public S reduce(T t) {
		return field.reduceInteger(t);
	}

	@Override
	public T upToUniformizerPower(T t) {
		return field.upToUniformizerPower(t);
	}

	@Override
	public Polynomial<S> reducePolynomial(Polynomial<T> t) {
		PolynomialRing<S> ring = AbstractPolynomialRing.getPolynomialRing(reduction(), t.numberOfVariables(),
				t.getPolynomialRing().getComparator());
		return ring.getEmbedding(t, new MathMap<>() {
			@Override
			public S evaluate(T t) {
				return reduce(t);
			}
		});
	}

	@Override
	public UnivariatePolynomial<S> reduceUnivariatePolynomial(UnivariatePolynomial<T> t) {
		return reduction().getUnivariatePolynomialRing().toUnivariate(reducePolynomial(t));
	}

	@Override
	public T lift(S s) {
		return field.liftToInteger(s);
	}

	@Override
	public Polynomial<T> liftPolynomial(Polynomial<S> s) {
		PolynomialRing<T> ring = AbstractPolynomialRing.getPolynomialRing(this, s.numberOfVariables(),
				s.getPolynomialRing().getComparator());
		return ring.getEmbedding(s, new MathMap<>() {
			@Override
			public T evaluate(S t) {
				return lift(t);
			}
		});
	}

	@Override
	public UnivariatePolynomial<T> liftUnivariatePolynomial(Polynomial<S> t) {
		return getUnivariatePolynomialRing().toUnivariate(liftPolynomial(t));
	}

	@Override
	public T round(T s, int accuracy) {
		return field.round(s, accuracy);
		// return getIdeal(Collections.singletonList(power(uniformizer(),
		// accuracy))).residue(s);
	}

	@Override
	public Polynomial<T> roundPolynomial(Polynomial<T> s, int accuracy) {
		PolynomialRing<T> ring = AbstractPolynomialRing.getPolynomialRing(this, s.numberOfVariables(),
				s.getPolynomialRing().getComparator());
		return ring.getEmbedding(s, new MathMap<>() {
			@Override
			public T evaluate(T t) {
				return round(t, accuracy);
			}
		});
	}

	@Override
	public UnivariatePolynomial<T> roundUnivariatePolynomial(Polynomial<T> t, int accuracy) {
		return getUnivariatePolynomialRing().toUnivariate(roundPolynomial(t, accuracy));
	}

	public boolean hasGoodReduction(UnivariatePolynomial<T> t) {
		if (!isIntegral(t)) {
			return false;
		}
		FactorizationResult<Polynomial<S>, S> reducedSquareFreeFactors = reduction().getUnivariatePolynomialRing()
				.squareFreeFactorization(reduceUnivariatePolynomial(t));
		return reducedSquareFreeFactors.squareFree();
	}

	public boolean hasIrreducibleGoodReduction(UnivariatePolynomial<T> t) {
		if (!isIntegral(t)) {
			return false;
		}
		return reduction().isIrreducible(reduceUnivariatePolynomial(t));
	}

	public boolean isEisenstein(UnivariatePolynomial<T> t) {
		if (t.degree() < 0) {
			return false;
		}
		if (valuation(t.leadingCoefficient()).value() != 0) {
			return false;
		}
		if (!valuation(t.univariateCoefficient(0)).equals(new Value(1))) {
			return false;
		}
		for (int i = 1; i < t.degree(); i++) {
			if (valuation(t.univariateCoefficient(i)).compareTo(new Value(0)) <= 0) {
				return false;
			}
		}
		return true;
	}

	public boolean isIntegral(UnivariatePolynomial<T> t) {
		if (valuation(t.leadingCoefficient()).value() != 0) {
			return false;
		}
		for (int i = 0; i < t.degree(); i++) {
			if (t.univariateCoefficient(i).equals(zero())) {
				continue;
			}
			if (valuation(t.univariateCoefficient(i)).compareTo(new Value(0)) < 0) {
				return false;
			}
		}
		return true;
	}

	public UnivariatePolynomial<T> integralMinimalPolynomial(UnivariatePolynomial<T> minimalPolynomial) {
		return integralPolynomial(minimalPolynomial).getFirst();
	}

	public UnivariatePolynomial<T> integralMinimalPolynomial(UnivariatePolynomial<T> minimalPolynomial, T uniformizer) {
		return integralPolynomial(minimalPolynomial, uniformizer).getFirst();
	}

	public Pair<UnivariatePolynomial<T>, UnivariatePolynomial<T>> integralPolynomial(UnivariatePolynomial<T> t) {
		return integralPolynomial(t, uniformizer());
	}

	public Pair<UnivariatePolynomial<T>, UnivariatePolynomial<T>> integralPolynomial(UnivariatePolynomial<T> t,
			T uniformizer) {
		UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
		int uniformizerValue = valuation(uniformizer).value();
		int lcValuation = valuation(t.leadingCoefficient()).value();
		int degree = t.degree();
		Value power = Value.INFINITY;
		for (int i = 0; i < degree; i++) {
			if (t.univariateCoefficient(i).equals(zero())) {
				continue;
			}
			int lowestPower = MiscAlgorithms.DivRoundDown(valuation(t.univariateCoefficient(i)).value() - lcValuation,
					uniformizerValue);
			power = power.min(new Value(MiscAlgorithms.DivRoundDown(lowestPower, degree - i)));
		}
		if (power.isInfinite()) {
			power = Value.ZERO;
		}
		T factor = field.power(uniformizer, power.value());
		UnivariatePolynomial<T> undoSubstitution = polynomials.divideScalar(polynomials.getVar(), factor);
		if (lcValuation + power.value() * t.degree() < 0) {
			t = polynomials.substitute(t,
					Collections.singletonList(polynomials.multiply(factor, polynomials.getVar())));
			t = polynomials.multiply(
					field.power(uniformizer,
							MiscAlgorithms.DivRoundDown(-lcValuation - power.value() * t.degree(), uniformizerValue)),
					t);
		} else {
			t = polynomials.multiply(
					field.power(uniformizer,
							MiscAlgorithms.DivRoundDown(-lcValuation - power.value() * t.degree(), uniformizerValue)),
					t);
			t = polynomials.substitute(t,
					Collections.singletonList(polynomials.multiply(factor, polynomials.getVar())));
		}
		return new Pair<>(polynomials.normalize(t), undoSubstitution);
	}

	@Override
	public T henselLift(UnivariatePolynomial<T> f, S aReduced) {
		return henselLift(f, aReduced, field.getAccuracy());
	}

	@Override
	public T henselLift(UnivariatePolynomial<T> f, S aReduced, int accuracy) {
		UnivariatePolynomial<S> fReduced = reduceUnivariatePolynomial(f);
		UnivariatePolynomialRing<S> reducedRing = reduction().getUnivariatePolynomialRing();
		if (fReduced.degree() != f.degree()) {
			throw new ArithmeticException("Leading coefficient divisible by p!");
		}
		if (!reducedRing.evaluate(fReduced, Collections.singletonList(aReduced)).equals(reduction().zero())) {
			throw new ArithmeticException("reduced polynomial does not have a zero at a!");
		}
		if (reducedRing.evaluate(reducedRing.derivative(fReduced), Collections.singletonList(aReduced))
				.equals(reduction().zero())) {
			throw new ArithmeticException("reduced derivative does have a zero at a!");
		}
		T a = lift(aReduced);
		return henselLiftWithInitialLift(f, a, accuracy);
	}

	@Override
	public T henselLiftWithInitialLift(UnivariatePolynomial<T> f, T initialLift) {
		return henselLiftWithInitialLift(f, initialLift, field.getAccuracy());
	}

	@SuppressWarnings("unchecked")
	@Override
	public T henselLiftWithInitialLift(fields.interfaces.UnivariatePolynomial<T> f, T initialLift, int accuracy) {
		UnivariatePolynomialRing<T> ring = getUnivariatePolynomialRing();
		UnivariatePolynomial<T> derivative = ring.derivative(f);
		Value initialAccuracy = valuation(ring.evaluate(f, initialLift));
		if (initialAccuracy.isInfinite()) {
			return initialLift;
		}
		int achievedAccuracy = initialAccuracy.value();
		if (achievedAccuracy <= 0) {
			throw new ArithmeticException("Initial lift not zero");
		}
		T a = initialLift;
		while (achievedAccuracy < accuracy) {
			a = subtract(a, divide(ring.evaluate(f, a), ring.evaluate(derivative, a)));
			achievedAccuracy *= 2;
		}
		return a;
	}

	@Override
	public UnivariatePolynomial<T> henselLiftFactor(UnivariatePolynomial<T> f, UnivariatePolynomial<S> gReduced) {
		return henselLiftFactor(f, gReduced, field.getAccuracy());
	}

	@Override
	public UnivariatePolynomial<T> henselLiftFactor(UnivariatePolynomial<T> f, UnivariatePolynomial<S> gReduced,
			int accuracy) {
		return extendedHenselLiftFactor(f, gReduced, accuracy).getLiftOfFactor();
	}

	@Override
	public HenselLiftFactorResult<T> extendedHenselLiftFactor(UnivariatePolynomial<T> f,
			UnivariatePolynomial<S> gReduced) {
		return extendedHenselLiftFactor(f, gReduced, field.getAccuracy());
	}

	@Override
	public HenselLiftFactorResult<T> extendedHenselLiftFactor(UnivariatePolynomial<T> f,
			UnivariatePolynomial<S> gReduced, int accuracy) {
		UnivariatePolynomialRing<T> ring = getUnivariatePolynomialRing();
		Polynomial<S> fReduced = reduceUnivariatePolynomial(f);
		if (fReduced.degree() != f.degree()) {
			throw new ArithmeticException("Leading coefficient divisible by p!");
		}
		if (!gReduced.leadingCoefficient().equals(reduction().one())) {
			throw new ArithmeticException("g not normalized");
		}
		if (!f.leadingCoefficient().equals(one())) {
			throw new ArithmeticException("f not normalized");
		}
		UnivariatePolynomialRing<S> reducedRing = reduction().getUnivariatePolynomialRing();
		// fReduced = reduceUnivariatePolynomial(f);
		Polynomial<S> hReduced = reducedRing.divideChecked(fReduced, gReduced);
		ExtendedEuclideanResult<Polynomial<S>> ee = reducedRing.extendedEuclidean(gReduced, hReduced);
		if (ee.getGcd().degree() != 0) {
			throw new ArithmeticException("Not coprime! f: " + f + " gReduced: " + gReduced + " gcd: " + ee.getGcd());
		}
		UnivariatePolynomial<T> g = liftUnivariatePolynomial(gReduced);
		Polynomial<T> h = liftPolynomial(hReduced);
		S gcd = reduction().inverse(ee.getGcd().leadingCoefficient());
		Polynomial<S> aReduced = reducedRing.toUnivariate(reducedRing.multiply(gcd, ee.getCoeff1()));
		Polynomial<S> bReduced = reducedRing.multiply(gcd, ee.getCoeff2());
		Polynomial<T> a = liftPolynomial(aReduced);
		Polynomial<T> b = liftPolynomial(bReduced);
		int achievedAccuracy = 1;
		while (achievedAccuracy < accuracy) {
			achievedAccuracy *= 2;
			Polynomial<T> q = ring.subtract(f, ring.multiply(g, h));
			g = ring.add(g, ring.multiply(b, q));
			h = ring.add(h, ring.multiply(a, q));
			g = roundUnivariatePolynomial(g, achievedAccuracy);
			h = roundPolynomial(h, achievedAccuracy);
			while (g.degree() != gReduced.degree() || !g.leadingCoefficient().equals(one())) {
				Polynomial<T> multiplier = null;
				if (g.degree() != gReduced.degree()) {
					multiplier = ring.getEmbedding(
							divide(g.leadingCoefficient(), g.univariateCoefficient(gReduced.degree())),
							g.degree() - gReduced.degree());
				} else {
					multiplier = ring.getEmbedding(subtract(one(), inverse(g.leadingCoefficient())));
				}
				g = ring.multiply(g, ring.subtract(ring.one(), multiplier));
				h = ring.multiply(h, ring.add(ring.one(), multiplier));
				g = roundUnivariatePolynomial(g, achievedAccuracy);
				h = roundPolynomial(h, achievedAccuracy);
			}
			Polynomial<T> r = ring.add(ring.multiply(a, g), ring.multiply(b, h), ring.negative(ring.one()));
			r = roundPolynomial(r, achievedAccuracy);
			a = ring.multiply(a, ring.subtract(ring.one(), r));
			b = ring.multiply(b, ring.subtract(ring.one(), r));
			a = roundPolynomial(a, achievedAccuracy);
			b = roundPolynomial(b, achievedAccuracy);
			while (a.degree() >= h.degree()) {
				Polynomial<T> multiplier = ring.getEmbedding(a.leadingCoefficient(), a.degree() - h.degree());
				a = ring.subtract(a, ring.multiply(multiplier, h));
				b = ring.add(b, ring.multiply(multiplier, g));
				a = roundPolynomial(a, achievedAccuracy);
				b = roundPolynomial(b, achievedAccuracy);
			}
		}
		return new HenselLiftFactorResult<>(g, ring.toUnivariate(h), ring.toUnivariate(a), ring.toUnivariate(b));
	}

	@Override
	public FactorizationResult<Polynomial<T>, T> henselLiftFactorization(UnivariatePolynomial<T> f, int accuracy) {
		UnivariatePolynomial<S> fReduced = reduceUnivariatePolynomial(f);
		if (fReduced.degree() != f.degree()) {
			throw new ArithmeticException("Leading coefficient divisible by p!");
		}
		UnivariatePolynomialRing<T> polynomialRing = getUnivariatePolynomialRing();
		UnivariatePolynomialRing<S> reducedPolynomialRing = reduction().getUnivariatePolynomialRing();
		FactorizationResult<Polynomial<S>, S> reducedFactorization = reduction().factorization(fReduced);
		List<Polynomial<T>> factors = new ArrayList<>();
		T unit = lift(reducedFactorization.getUnit());
		Polynomial<T> product = polynomialRing.getEmbedding(unit);
		List<Polynomial<S>> cofactors = new ArrayList<>();
		for (Polynomial<S> factor : reducedFactorization.primeFactors()) {
			if (reducedFactorization.multiplicity(factor) != 1) {
				throw new ArithmeticException("Reduction not square free");
			}
			Polynomial<T> liftedFactor = liftPolynomial(factor);
			factors.add(liftedFactor);
			cofactors.add(reducedPolynomialRing.divideChecked(fReduced, factor));
			product = polynomialRing.multiply(product, liftedFactor);
		}
		List<Polynomial<T>> initialFactors = new ArrayList<>();
		initialFactors.addAll(factors);
		ExtendedEuclideanListResult<Polynomial<S>> gcd = reducedPolynomialRing.extendedEuclidean(cofactors);
		List<Polynomial<T>> delta = new ArrayList<>();
		for (Polynomial<S> coefficient : gcd.getCoeffs()) {
			delta.add(liftPolynomial(coefficient));
		}
		for (int k = 1; k < accuracy; k++) {
			Polynomial<T> error = roundPolynomial(polynomialRing.subtract(f, product), k + 1);
			product = polynomialRing.getEmbedding(unit);
			for (int i = 0; i < factors.size(); i++) {
				Polynomial<T> factor = factors.get(i);
				Polynomial<T> adjustment = roundPolynomial(
						polynomialRing.remainder(polynomialRing.multiply(error, delta.get(i)), initialFactors.get(i)),
						k + 1);
				Polynomial<T> nextFactor = polynomialRing.add(factor, adjustment);
				product = roundPolynomial(polynomialRing.multiply(product, nextFactor), k + 2);
				factors.set(i, nextFactor);
			}
		}
		SortedMap<Polynomial<T>, Integer> factorization = new TreeMap<>();
		for (Polynomial<T> factor : factors) {
			factorization.put(factor, 1);
		}
		return new FactorizationResult<>(unit, factorization);
	}

	@Override
	public String toString() {
		return asString;
	}

	@Override
	public T characteristicRoot(T t, int power) {
		return field.characteristicRoot(t, power);
	}

	@Override
	public boolean hasCharacteristicRoot(T t, int power) {
		return field.hasCharacteristicRoot(t, power);
	}

	@Override
	public boolean isIrreducible(UnivariatePolynomial<T> t) {
		return theMontesAlgorithm(t).getTypes().size() == 1;
//		t = integralMinimalPolynomial(t);
//		if (isEisenstein(t)) {
//			return true;
//		}
//		UnivariatePolynomial<S> reduced = reduceUnivariatePolynomial(t);
//		if (reduced.degree() == t.degree() && reduction().isIrreducible(reduced)) {
//			return true;
//		}
//		UnivariatePolynomialRing<S> reductionPolynomials = reduction().getUnivariatePolynomialRing();
//		Map<Polynomial<S>, Integer> factors = reductionPolynomials.squareFreeFactorization(reduced);
//		if (factors.size() != 1) {
//			return false;
//		}
//		Polynomial<S> uniqueFactor = factors.keySet().iterator().next();
//		if (!reduction().isIrreducible(reductionPolynomials.toUnivariate(uniqueFactor))) {
//			return false;
//		}
//		UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
//		Polynomial<T> lifted = liftPolynomial(uniqueFactor);
//		CoordinateRing<T> cr = new CoordinateRing<>(polynomials, polynomials.getIdeal(Collections.singletonList(t)));
//		int ramificationIndex = t.degree() / lifted.degree();
//		CoordinateRingElement<T> liftedCr = cr.getEmbedding(lifted);
//		CoordinateRingElement<T> uniformizerInverted = cr.getEmbedding(field.inverse(field.uniformizer()));
//		int val = 1;
//		while (true) {
//			CoordinateRingElement<T> probe = cr.multiply(cr.power(liftedCr, ramificationIndex),
//					cr.power(uniformizerInverted, val));
//			if (probe.equals(cr.zero())) {
//				return false;
//			}
//			probe = cr.power(probe, 8);
//			Value lowestValuation = new Value(0);
//			for (int j = 0; j <= probe.getElement().degree(); j++) {
//				T c = polynomials.toUnivariate(probe.getElement()).univariateCoefficient(j);
//				Value valuation = field.valuation(c);
//				if (lowestValuation.compareTo(valuation) > 0) {
//					lowestValuation = valuation;
//				}
//			}
//			if (lowestValuation.compareTo(new Value(-val)) < 0) {
//				return false;
//			}
//			probe = cr.power(probe, field.getAccuracy());
//			if (!probe.equals(cr.zero())) {
//				return true;
//			}
//			val++;
//		}
	}

	private class Type<R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>>
			implements OkutsuType<T, S, R, RE, RFE>, Cloneable {
		private UnivariatePolynomial<T> polynomial;
		private int level;
		private Type<R, RE, RFE> prevLevel;
		private final int order;
		private final int slopes;
		private UnivariatePolynomial<T> representative;
		private UnivariatePolynomial<RE> reduced;
		private int value;
		private Fraction nu;
		private Fraction precision;
		private Fraction lambda;
		private Fraction w;
		private int ramificationIndex;
		private int residueDegree;
		private int h;
		private int e;
		private int f;
		private int l;
		private int lPrime;
		private FieldEmbedding<R, RE, RFE> reductionFieldExtension;
		private Extension<S, R, RE, RFE> reductionAsExtension;
		private int[] logUniformizer;
		private int[] capitalPhi;
		private int[] gamma;
		private Map<Integer, UnivariatePolynomial<T>> divisorPolynomials;
		private OkutsuType<T, S, R, RE, RFE> doublePrecision;

		// Constructor for level 0
		private Type(UnivariatePolynomial<T> polynomial, UnivariatePolynomial<T> representative, int order,
				Extension<S, R, RE, RFE> reductionAsExtension) {
			this.level = 0;
			this.polynomial = polynomial;
			this.ramificationIndex = 1;
			this.lambda = Rationals.q().zero();
			this.nu = Rationals.q().zero();
			this.precision = Rationals.q().zero();
			this.w = Rationals.q().zero();
			this.h = 0;
			this.e = 1;
			this.representative = representative;
			this.order = order;
			this.slopes = 1;
			UnivariatePolynomial<S> reduced = reduceUnivariatePolynomial(representative);
			RFE reduction = reductionAsExtension.extension();
			UnivariatePolynomialRing<RE> reductionPolynomials = reduction.getUnivariatePolynomialRing();
			UnivariatePolynomial<RE> reducedEmbedded = reductionPolynomials.getEmbedding(reduced,
					reductionAsExtension.embeddingMap());
			this.reduced = reducedEmbedded;
			this.residueDegree = reduced.degree();
			this.f = this.residueDegree;
			this.l = 0;
			this.lPrime = 1;
			this.value = 0;
			this.reductionFieldExtension = reduction.getEmbeddedExtension(reducedEmbedded);
			this.reductionAsExtension = reductionAsExtension;
			this.logUniformizer = new int[] { 1 };
			this.capitalPhi = new int[] { 0 };
			this.gamma = new int[] { 0 };
		}

		// Constructor for divisors
		private Type(Type<R, RE, RFE> prevLevel, UnivariatePolynomial<T> representative) {
			this.prevLevel = prevLevel;
			this.representative = representative;
			this.slopes = 1;
			this.order = 1;
			computeDerivedValues();
		}

		// Constructor for non root nodes
		private Type(Type<R, RE, RFE> prevLevel, UnivariatePolynomial<T> representative,
				UnivariatePolynomial<RE> reduced, Fraction lambda, int order, int slopes) {
			this.prevLevel = prevLevel;
			this.representative = representative;
			this.reduced = reduced;
			this.lambda = lambda;
			this.order = order;
			this.slopes = slopes;
			computeDerivedValues();
		}

		// Copy Constructor
		private Type(Type<R, RE, RFE> copy) {
			this.level = copy.level;
			this.prevLevel = copy.prevLevel == null ? null : copy.prevLevel.clone();
			this.polynomial = copy.polynomial;
			this.ramificationIndex = copy.ramificationIndex;
			this.lambda = copy.lambda;
			this.nu = copy.nu;
			this.precision = copy.precision;
			this.w = copy.w;
			this.h = copy.h;
			this.e = copy.e;
			this.representative = copy.representative;
			this.order = copy.order;
			this.slopes = copy.slopes;
			this.reduced = copy.reduced;
			this.residueDegree = copy.residueDegree;
			this.f = copy.f;
			this.l = copy.l;
			this.lPrime = copy.lPrime;
			this.value = copy.value;
			this.reductionFieldExtension = copy.reductionFieldExtension;
			this.reductionAsExtension = copy.reductionAsExtension;
			this.logUniformizer = copy.logUniformizer;
			this.capitalPhi = copy.capitalPhi;
			this.gamma = copy.gamma;
		}

		private void computeDerivedValues() {
			Integers z = Integers.z();
			Rationals q = Rationals.q();
			this.level = prevLevel.level + 1;
			this.polynomial = prevLevel.polynomial;
			if (lambda == null) {
				this.value = -1;
				this.nu = null;
				this.w = null;
				this.precision = null;
				this.e = 1;
				this.f = 1;
				this.l = 0;
				this.lPrime = 1;
				this.reduced = prevLevel.reductionFieldExtension.getField().getUnivariatePolynomialRing().getVar();
			} else {
				this.h = this.lambda.getNumerator().getValue().intValueExact();
				this.e = this.lambda.getDenominator().getValue().intValueExact();
				this.f = reduced.degree();
				ExtendedEuclideanResult<IntE> ee = z.extendedEuclidean(z.getInteger(h), z.getInteger(e));
				this.l = ee.getCoeff1().getValue().intValueExact();
				this.lPrime = ee.getCoeff2().getValue().intValueExact();
				if (this.l < 0) {
					this.l += e;
					this.lPrime -= h;
				}
			}
			this.ramificationIndex = prevLevel.ramificationIndex * e;
			this.residueDegree = prevLevel.residueDegree * f;
			if (lambda != null) {
				this.nu = q.getFraction(z.getInteger(h), z.getInteger(ramificationIndex));
				this.precision = q.add(this.nu, prevLevel.precision);
				this.w = q.multiply(q.getInteger(e * f), q.add(prevLevel.w, nu));
				this.value = q.multiply(this.w, q.getInteger(ramificationIndex)).getNumerator().getValue()
						.intValueExact();
				this.capitalPhi = new int[level + 1];
				for (int i = 0; i < this.capitalPhi.length - 1; i++) {
					this.capitalPhi[i] = -prevLevel.logUniformizer[i] * prevLevel.value;
				}
				this.capitalPhi[level] = 1;
				this.gamma = new int[level + 1];
				for (int i = 0; i < this.gamma.length - 1; i++) {
					this.gamma[i] = this.capitalPhi[i] * e - prevLevel.logUniformizer[i] * h;
				}
				this.gamma[level] = this.capitalPhi[level] * e;
				this.logUniformizer = new int[level + 1];
				for (int i = 0; i < this.logUniformizer.length - 1; i++) {
					this.logUniformizer[i] = this.capitalPhi[i] * l + prevLevel.logUniformizer[i] * lPrime;
				}
				this.logUniformizer[level] = this.capitalPhi[level] * l;
			}
			this.reductionFieldExtension = prevLevel.reductionFieldExtension.getField()
					.getEmbeddedExtension(this.reduced);
			this.reductionAsExtension = new Extension<>(reductionFieldExtension.getField(),
					prevLevel.reductionAsExtension.domain(),
					new ConCatMap<>(prevLevel.reductionAsExtension.embeddingMap(),
							reductionFieldExtension.getEmbeddingMap()),
					new ConCatMap<>(reductionFieldExtension.asVectorMap(),
							new FlattenProductVectorMap<>(prevLevel.reductionAsExtension.asVectorMap())));
		}

		private int[] phiToGamma(int[] c) {
			int[] vector = Arrays.copyOf(c, c.length);
			int[] result = new int[level];
			for (Type<R, RE, RFE> it = this; it.level > 0; it = it.prevLevel) {
				result[it.level - 1] = vector[it.level] / it.e;
				for (int i = 0; i < it.gamma.length; i++) {
					vector[i] -= result[it.level - 1] * it.gamma[i];
				}
			}
			return result;
		}

		@Override
		public Type<R, RE, RFE> clone() {
			return new Type<>(this);
		}

		private Pair<Integer, UnivariatePolynomial<T>> computeDivisorPolynomial(int m) {
			UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
			if (level == 0) {
				int base = representative.degree();
				int c = m % base;
				return new Pair<>((m - c) / base, polynomials.getVarPower(c));
			}
			Pair<Integer, UnivariatePolynomial<T>> prev = prevLevel.computeDivisorPolynomial(m);
			m = prev.getFirst();
			UnivariatePolynomial<T> divisorPolynomial = prev.getSecond();
			int base = representative().degree() / phi().degree();
			int c = m % base;
			divisorPolynomial = polynomials.multiply(divisorPolynomial, polynomials.power(phi(), c));
			return new Pair<>((m - c) / base, divisorPolynomial);
		}

		@Override
		public UnivariatePolynomial<T> divisorPolynomial(int m) {
			if (!leaf()) {
				throw new ArithmeticException("Divisor Polynomial only defined at leaf!");
			}
			if (m >= representative().degree()) {
				throw new ArithmeticException("m too large");
			}
			if (divisorPolynomials == null) {
				divisorPolynomials = new TreeMap<>();
			} else if (divisorPolynomials.containsKey(m)) {
				return divisorPolynomials.get(m);
			}
			Pair<Integer, UnivariatePolynomial<T>> result = prevLevel.computeDivisorPolynomial(m);
			if (result.getFirst() != 0) {
				throw new ArithmeticException("m cannot be expressed");
			}
			divisorPolynomials.put(m, result.getSecond());
			return result.getSecond();
		}

		@Override
		public String toString() {
			if (level == 0) {
				return "L0, Rep: " + representative;
			}
			return "L" + level + " (" + prevLevel.representative + ", " + lambda + ", " + reduced + "), Rep: "
					+ representative;
		}

		@Override
		public UnivariatePolynomial<T> getPolynomial() {
			return polynomial;
		}

		@Override
		public UnivariatePolynomial<T> representative() {
			return representative;
		}

		@Override
		public UnivariatePolynomial<T> phi() {
			return prevLevel.representative;
		}

		@Override
		public Fraction lambda() {
			return lambda;
		}

		@Override
		public UnivariatePolynomial<RE> psi() {
			return reduced;
		}

		@Override
		public int level() {
			return level;
		}

		@Override
		public OkutsuType<T, S, R, RE, RFE> previousLevel() {
			return prevLevel;
		}

		@Override
		public boolean complete() {
			return order == 1 && slopes == 1;
		}

		@Override
		public boolean leaf() {
			return (level > 0 && prevLevel.complete()) || lambda == null;
		}

		@Override
		public Value value() {
			if (lambda == null) {
				return Value.INFINITY;
			}
			if (level == 0) {
				return Value.ZERO;
			}
			return new Value(prevLevel.value);
		}

		@Override
		public Value quality() {
			if (lambda == null) {
				return Value.INFINITY;
			}
			return new Value(value);
		}

		@Override
		public Value precision() {
			if (lambda == null) {
				return Value.INFINITY;
			}
			return new Value(MiscAlgorithms.DivRoundUp(precision.getNumerator().getValue().intValueExact(),
					precision.getDenominator().getValue().intValueExact()));
		}

		@Override
		public int exponent() {
			if (lambda == null) {
				return -1;
			}
			if (!leaf()) {
				throw new ArithmeticException("not a leaf");
			}
			Rationals q = Rationals.q();
			Fraction mu = q.subtract(q.getFraction(new IntE(value), new IntE(ramificationIndex)), prevLevel.precision);
			return MiscAlgorithms.DivRoundDown(mu.getNumerator().getValue().intValueExact(),
					mu.getDenominator().getValue().intValueExact());
		}

		@Override
		public int ramificationIndex() {
			return ramificationIndex;
		}

		@Override
		public int residueDegree() {
			return residueDegree;
		}

		@Override
		public int[] uniformizer() {
			return logUniformizer;
		}

		@Override
		public Extension<S, R, RE, RFE> reduction() {
			return reductionAsExtension;
		}

		@Override
		public FieldEmbedding<R, RE, RFE> embeddedReductionExtension() {
			return reductionFieldExtension;
		}

		@Override
		public Type<R, RE, RFE> updateLeaf(UnivariatePolynomial<T> representative) {
			if (!complete()) {
				throw new ArithmeticException("Not a complete type!");
			}
			this.doublePrecision = null;
			if (leaf()) {
				Type<R, RE, RFE> result = prevLevel.updateLeaf(representative);
				if (result.precision().compareTo(precision()) < 0) {
					return this;
//					throw new ArithmeticException("Not a precision upgrade!");
				}
				this.prevLevel = result.prevLevel;
				this.representative = result.representative;
				this.reduced = result.reduced;
				this.lambda = result.lambda;
				computeDerivedValues();
				return result;
			}
			if (representative.degree() != this.representative.degree()) {
				throw new ArithmeticException("different degree in representative!");
			}
			Type<R, RE, RFE> copy = clone();
			copy.representative = representative;
			if (field.getUnivariatePolynomialRing().isDivisible(getPolynomial(), representative)) {
				return new Type<>(copy, representative);
			}
			NewtonPolygon polygon = copy.newtonPolygon(getPolynomial(), 1);
			if (polygon.slopes.size() != 1) {
				throw new ArithmeticException("Expected exactly one slope");
			}
			if (polygon.coefficients.size() != 2) {
				throw new ArithmeticException("Expected exactly two points in Newton polygon");
			}
			Fraction lambda = Rationals.q().negative(polygon.slopes.get(0));
			UnivariatePolynomial<RE> reduced = reductionFieldExtension.getField().getUnivariatePolynomialRing()
					.normalize(copy.reduceNextLevel(getPolynomial(), lambda, polygon));
			UnivariatePolynomial<T> nextRepresentative = copy.representativeNextLevel(reduced, lambda);
			return new Type<>(copy, nextRepresentative, reduced, lambda, 1, 1);
		}

		@Override
		public UnivariatePolynomial<T> liftAsPolynomial(UnivariatePolynomial<RE> reduced, int b) {
			UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
			if (level == 0) {
				UnivariatePolynomialRing<S> reductionPolynomials = LocalRingImplementation.this.reduction()
						.getUnivariatePolynomialRing();
				UnivariatePolynomial<S> asReduced = reductionPolynomials.getEmbedding(reduced,
						reductionAsExtension.retractionMap());
				UnivariatePolynomial<T> lifted = liftUnivariatePolynomial(asReduced);
				return polynomials.multiply(field.power(field.uniformizer(), b), lifted);
			}
			if (lambda != null && b < value) {
				int add = MiscAlgorithms.DivRoundUp(value - b, ramificationIndex);
				RE multiply = reduce(polynomials.getEmbedding(field.power(field.uniformizer(), add)));
				RE asElement = reductionFieldExtension.fromPolynomial(reduced);
				RE toLift = reductionFieldExtension.getField().multiply(multiply, asElement);
				return polynomials.multiply(field.power(field.uniformizer(), -add),
						lift(toLift, b + add * ramificationIndex));
			}
			return prevLevel.liftNextLevel(reduced, lambda, b);
		}

		@Override
		public UnivariatePolynomial<T> lift(RE reduced, int b) {
			return liftAsPolynomial(reductionFieldExtension.asPolynomial(reduced), b);
		}

		private UnivariatePolynomial<T> liftNextLevel(UnivariatePolynomial<RE> reduced, Fraction lambda, int b) {
			RFE field = reductionFieldExtension.getField();
			if (reduced.equals(field.getUnivariatePolynomialRing().zero())) {
				return getUnivariatePolynomialRing().zero();
			}
			Rationals q = Rationals.q();
			IntE h = lambda == null ? null : lambda.getNumerator();
			IntE e = lambda == null ? Integers.z().one() : lambda.getDenominator();
			int s = h == null ? 0
					: h.getValue().modInverse(e.getValue()).multiply(BigInteger.valueOf(b)).mod(e.getValue())
							.intValueExact();
			int eMinusOne = this.e;
			int lMinusOne = this.l;
			// int value = q.multiply(q.getInteger(e * this.f * this.ramificationIndex),
			// q.add(this.nu, this.w)).getNumerator().getValue().intValueExact();

			int k = reduced.order().value();
			UnivariatePolynomialRing<T> polynomials = getUnivariatePolynomialRing();
			UnivariatePolynomial<T> liftedPower = polynomials
					.toUnivariate(polynomials.power(representative, s + k * e.getValue().intValueExact()));
			UnivariatePolynomial<T> lifted = polynomials.zero();
			UnivariatePolynomial<T> base = polynomials
					.toUnivariate(polynomials.power(representative, e.getValue().intValueExact()));
			UnivariatePolynomial<T> basePower = polynomials.one();
			for (int j = 0; j <= reduced.degree() - k; j++) {
				int sj = s + (j + k) * e.getValue().intValueExact();
				Fraction bjFraction = q.subtract(q.divide(q.getInteger(b), q.getInteger(e)),
						lambda == null ? q.zero() : q.multiply(sj, q.add(q.getInteger(value), lambda)));
				bjFraction.canonicalize();
				if (!bjFraction.getDenominator().equals(Integers.z().one())) {
					throw new ArithmeticException("bj computation wrong!");
				}
				int bj = bjFraction.getNumerator().getValue().intValueExact();
				int adjustedSj = (lMinusOne * bj) % eMinusOne;
				int power = (lMinusOne * bj - adjustedSj) / eMinusOne;
				RE element = field.multiply(reduced.univariateCoefficient(j + k),
						field.power(reductionFieldExtension.getGenerator(), power));
				UnivariatePolynomial<T> lift = lift(element, bj);
				lifted = polynomials.add(lifted, polynomials.multiply(lift, basePower));
				basePower = polynomials.multiply(basePower, base);
			}
			return polynomials.multiply(lifted, liftedPower);
		}

		private UnivariatePolynomial<T> representativeNextLevel(UnivariatePolynomial<RE> reduced, Fraction lambda) {
			RFE field = reductionFieldExtension.getField();
			UnivariatePolynomialRing<RE> reductionPolynomials = field.getUnivariatePolynomialRing();
			UnivariatePolynomial<RE> adjustedReduced = reductionPolynomials.toUnivariate(
					reductionPolynomials.subtract(reduced, reductionPolynomials.getVarPower(reduced.degree())));
			Rationals q = Rationals.q();
			int e = lambda.getDenominator().getValue().intValueExact();
			int f = reduced.degree();
			// int requiredValue = e * f * value;
			Fraction nu = q.getFraction(lambda.getNumerator(), Integers.z().getInteger(e * this.ramificationIndex));
			int value = q.multiply(q.getInteger(e * f * e * this.ramificationIndex), q.add(nu, this.w)).getNumerator()
					.getValue().intValueExact();
			UnivariatePolynomial<T> lifted = liftNextLevel(adjustedReduced, lambda, value);
			UnivariatePolynomialRing<T> polynomials = getUnivariatePolynomialRing();
			return polynomials.add(polynomials.power(representative, e * f), lifted);
		}

		private OkutsuType<T, S, R, RE, RFE> doublePrecision() {
			if (doublePrecision == null) {
				Type<R, RE, RFE> copy = clone();
				doublePrecision = singleFactorLifting(copy, 2 * precision().value());
			}
			return doublePrecision;
		}

		@Override
		public Value valuation(UnivariatePolynomial<T> t) {
			if (level == 0) {
				return LocalRingImplementation.this.valuationOfUnivariatePolynomial(t);
			}
			Value prevLevelValue = prevLevel.valuationNextLevel(t, lambda);
//			if (!prevLevelValue.isInfinite() && this.leaf() && !quality().isInfinite()
//					&& t.degree() >= representative().degree()) {
//				int levelZero = LocalRingImplementation.this.valuationOfUnivariatePolynomial(t).value();
//				if (levelZero >= 0) {
//					levelZero = 0;
//				}
//				if (new Value(prevLevelValue.value() - levelZero * ramificationIndex()).compareTo(quality()) >= 0) {
			if (prevLevelValue == null) {
				return doublePrecision().valuation(t);
			}
//			}
			return prevLevelValue;
		}

		private Value valuationNextLevel(UnivariatePolynomial<T> t, Fraction lambda) {
			if (t.equals(getUnivariatePolynomialRing().zero())) {
				return Value.INFINITY;
			}
			Rationals q = Rationals.q();
			NewtonPolygon polygon = newtonPolygon(t, -1);
			if (lambda == null) {
				return polygon.value(0);
			}
			Pair<Integer, Integer> criticalLine = polygon.criticalLine(lambda);
			if (complete() && criticalLine.getFirst() != criticalLine.getSecond()) {
				return null;
			}
			Fraction s = q.getInteger(criticalLine.getFirst());
			Fraction v = q.getInteger(polygon.value(criticalLine.getFirst()).value());
			Fraction h = q.add(q.multiply(s, lambda), v);
			IntE e = lambda.getDenominator();
			Fraction value = q.multiply(q.getInteger(e), h);
			value.canonicalize();
			if (!value.getDenominator().equals(Integers.z().one())) {
				throw new ArithmeticException("Value not computed correctly");
			}
			return new Value(value.getNumerator().getValue().intValueExact());
		}

		@Override
		public UnivariatePolynomial<RE> reduceAsPolynomial(UnivariatePolynomial<T> t) {
			return reduceAndData(t).reduction;
		}

		@Override
		public RE reduce(UnivariatePolynomial<T> t) {
			return reductionFieldExtension.fromPolynomial(reduceAsPolynomial(t));
		}

		private class ReductionValuationResult {
			private UnivariatePolynomial<RE> reduction;
			private int power;

			public ReductionValuationResult(UnivariatePolynomial<RE> reduction, int power) {
				this.reduction = reduction;
				this.power = power;
			}
		}

		private ReductionValuationResult reduceAndData(UnivariatePolynomial<T> t) {
			if (level == 0) {
				UnivariatePolynomialRing<RE> reductionPolynomials = reductionFieldExtension.getEmbeddedField()
						.getUnivariatePolynomialRing();
				if (t.equals(getUnivariatePolynomialRing().zero())) {
					return new ReductionValuationResult(reductionPolynomials.zero(), 0);
				}
				int value = valuationOfUnivariatePolynomial(t).value();
				t = field.getUnivariatePolynomialRing().divideScalar(t, field.power(field.uniformizer(), value));
				return new ReductionValuationResult(reductionPolynomials.getEmbedding(reduceUnivariatePolynomial(t),
						reductionAsExtension.embeddingMap()), 0);
			}
			NewtonPolygon polygon = prevLevel.newtonPolygon(t, -1);
			return prevLevel.reduceAndDataNextLevel(t, lambda, polygon);
		}

		private UnivariatePolynomial<RE> reduceNextLevel(UnivariatePolynomial<T> t, Fraction lambda,
				NewtonPolygon polygon) {
			return reduceAndDataNextLevel(t, lambda, polygon).reduction;
		}

		private ReductionValuationResult reduceAndDataNextLevel(UnivariatePolynomial<T> t, Fraction lambda,
				NewtonPolygon polygon) {
			UnivariatePolynomialRing<T> polynomials = getUnivariatePolynomialRing();
			RFE field = reductionFieldExtension.getField();
			RE generator = reductionFieldExtension.getGenerator();
			UnivariatePolynomialRing<RE> reductionPolynomials = field.getUnivariatePolynomialRing();
			if (t.equals(getUnivariatePolynomialRing().zero())) {
				return new ReductionValuationResult(reductionPolynomials.zero(), 0);
			}
			Rationals q = Rationals.q();
			List<RE> result = new ArrayList<>();
			Pair<Integer, Integer> criticalLine = polygon.criticalLine(lambda);
			int s = criticalLine.getFirst();
			int sPrime = criticalLine.getSecond();
			int e = lambda == null ? 1 : lambda.getDenominator().getValue().intValueExact();
			Fraction startValue = q.getInteger(polygon.value(s).value());
			for (int i = s; i <= sPrime; i += e) {
				Fraction maxValue = q.add(startValue, s == i ? q.zero() : q.multiply(lambda, q.getInteger(s - i)));
				if (polygon.value(i).isInfinite() || q.getInteger(polygon.value(i).value()).compareTo(maxValue) > 0) {
					result.add(reductionFieldExtension.getField().zero());
					continue;
				}
				ReductionValuationResult reduceAndData = reduceAndData(
						polynomials.toUnivariate(polygon.coefficient(i)));
				RE power = field.power(generator, reduceAndData.power);
				UnivariatePolynomial<RE> reduced = reductionPolynomials.getEmbedding(reduceAndData.reduction,
						reductionFieldExtension.getEmbeddingMap());
				@SuppressWarnings("unchecked")
				RE eval = reductionPolynomials.evaluate(reduced, generator);
				result.add(field.multiply(power, eval));
			}
			UnivariatePolynomial<RE> reduced = reductionPolynomials.getPolynomial(result);
			if (lambda == null) {
				return new ReductionValuationResult(reduced, 0);
			}
			Fraction h = q.add(q.multiply(s, lambda), startValue);
			int value = q.multiply(e, h).getNumerator().getValue().intValueExact();
			int l = Integers.z().extendedEuclidean(lambda.getNumerator(), lambda.getDenominator()).getCoeff1()
					.getValue().intValueExact();
			if (l < 0) {
				l += e;
			}
			return new ReductionValuationResult(reduced, (s - l * value) / e);
		}

		@Override
		public boolean divides(UnivariatePolynomial<T> t) {
			UnivariatePolynomial<RE> reduced = reduceAsPolynomial(t);
			UnivariatePolynomial<RE> reducedPolynomial = reduceAsPolynomial(polynomial);
			return reductionFieldExtension.getEmbeddedField().getUnivariatePolynomialRing()
					.isDivisible(reducedPolynomial, reduced);
		}

		private class NewtonPolygon {
			private List<Fraction> slopes;
			private List<Vector<IntE>> borderPoints;
			private List<Polynomial<T>> coefficients;
			private List<Value> values;

			public NewtonPolygon(UnivariatePolynomial<T> g, List<Polynomial<T>> coefficients, List<Value> values) {
				this.coefficients = coefficients;
				this.values = values;
				Integers z = Integers.z();
				Rationals q = Rationals.q();
				FreeModule<IntE> vectors = new FreeModule<>(z, 2);
				List<Vector<IntE>> points = new ArrayList<>();
				Vector<IntE> firstPoint = null;
				for (int i = 0; i < values.size(); i++) {
					if (values.get(i).isInfinite()) {
						continue;
					}
					List<IntE> coords = new ArrayList<>();
					coords.add(z.getInteger(i));
					coords.add(z.getInteger(values.get(i).value()));
					if (firstPoint == null) {
						firstPoint = new Vector<>(coords);
					}
					points.add(new Vector<>(coords));
				}
				final Vector<IntE> first = firstPoint;
				Collections.sort(points, new Comparator<Vector<IntE>>() {
					@Override
					public int compare(Vector<IntE> o1, Vector<IntE> o2) {
						Vector<IntE> diff1 = vectors.subtract(o1, first);
						Vector<IntE> diff2 = vectors.subtract(o2, first);
						boolean diff1zero = diff1.equals(vectors.zero());
						boolean diff2zero = diff2.equals(vectors.zero());
						if (diff1zero && diff2zero) {
							return 0;
						}
						if (diff1zero) {
							return -1;
						}
						if (diff2zero) {
							return 1;
						}
						Fraction slope1 = q.getFraction(diff1.get(2), diff1.get(1));
						Fraction slope2 = q.getFraction(diff2.get(2), diff2.get(1));
						int cmp = slope1.compareTo(slope2);
						if (cmp != 0) {
							return cmp;
						}
						return o2.get(1).compareTo(o1.get(1));
					}
				});
				slopes = new ArrayList<>();
				borderPoints = new ArrayList<>();
				Deque<Vector<IntE>> pointStack = new LinkedList<>();
				int s = -1;
				for (Vector<IntE> point : points) {
					int pointS = point.get(1).getValue().intValueExact();
					if (pointS < s) {
						continue;
					}
					s = pointS;
					if (pointStack.isEmpty()) {
						pointStack.push(point);
						continue;
					}
					Vector<IntE> top = pointStack.pop();
					while (!pointStack.isEmpty() && !counterclockwise(pointStack.peek(), top, point)) {
						top = pointStack.pop();
					}
					pointStack.push(top);
					pointStack.push(point);
				}
				Vector<IntE> prevPoint = null;
				while (!pointStack.isEmpty()) {
					Vector<IntE> point = pointStack.pollLast();
					this.borderPoints.add(point);
					if (prevPoint != null) {
						Fraction slope = q.getFraction(z.subtract(point.get(2), prevPoint.get(2)),
								z.subtract(point.get(1), prevPoint.get(1)));
						if (slope.compareTo(q.zero()) > 0) {
							break;
						}
						this.slopes.add(slope);
					}
					prevPoint = point;
				}
			}

			private boolean counterclockwise(Vector<IntE> point1, Vector<IntE> point2, Vector<IntE> point3) {
				Integers z = Integers.z();
				IntE diffx21 = z.subtract(point2.get(1), point1.get(1));
				IntE diffy31 = z.subtract(point3.get(2), point1.get(2));
				IntE diffx31 = z.subtract(point3.get(1), point1.get(1));
				IntE diffy21 = z.subtract(point2.get(2), point1.get(2));
				IntE prod1 = z.multiply(diffx21, diffy31);
				IntE prod2 = z.multiply(diffx31, diffy21);
				return z.subtract(prod1, prod2).compareTo(z.zero()) > 0;
			}

			@Override
			public String toString() {
				if (borderPoints == null || borderPoints.isEmpty()) {
					return "()";
				}
				StringBuilder b = new StringBuilder();
				b.append("(");
				b.append(borderPoints.get(0));
				for (int i = 0; i < slopes.size(); i++) {
					b.append(" ");
					b.append(slopes.get(i));
					b.append(" ");
					b.append(borderPoints.get(i + 1));
				}
				b.append(")");
				return b.toString();
			}

			public Pair<Integer, Integer> criticalLine(Fraction lambda) {
				if (lambda == null) {
					return new Pair<>(0, 0);
				}
				int index = Collections.binarySearch(slopes, Rationals.q().negative(lambda));
				if (index >= 0) {
					int s = borderPoints.get(index).get(1).getValue().intValueExact();
					int sPrime = borderPoints.get(index + 1).get(1).getValue().intValueExact();
					return new Pair<>(s, sPrime);
				}
				index = -(index + 1);
				int s = borderPoints.get(index).get(1).getValue().intValueExact();
				return new Pair<>(s, s);
			}

			public Value value(int index) {
				return values.get(index);
			}

			public Polynomial<T> coefficient(int index) {
				return coefficients.get(index);
			}
		}

		private NewtonPolygon newtonPolygon(UnivariatePolynomial<T> g, int order) {
			UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
			Polynomial<T> remaining = g;
			Polynomial<T> power = polynomials.one();
			List<Polynomial<T>> coefficients = new ArrayList<>();
			List<Value> values = new ArrayList<>();
			for (int i = 0; i <= order || (order == -1 && !remaining.equals(polynomials.zero())); i++) {
				QuotientAndRemainderResult<Polynomial<T>> qr = polynomials.quotientAndRemainder(remaining,
						representative);
				remaining = qr.getQuotient();
				Polynomial<T> coefficient = qr.getRemainder();
				Value value = valuation(polynomials.multiply(power, coefficient));
				coefficients.add(coefficient);
				values.add(value);
				power = polynomials.multiply(power, representative);
			}
			return new NewtonPolygon(g, coefficients, values);
		}
	}

	private <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> List<OkutsuType<T, S, R, RE, RFE>> theMontesAlgorithm(
			UnivariatePolynomial<T> t, Type<R, RE, RFE> prevLevel) {
		if (prevLevel.leaf()) {
			return Collections.singletonList(prevLevel.updateLeaf(prevLevel.representative()));
		}
		List<Type<R, RE, RFE>> workStack = new ArrayList<>();
		List<OkutsuType<T, S, R, RE, RFE>> result = new ArrayList<>();
		UnivariatePolynomialRing<T> polynomials = getUnivariatePolynomialRing();
		QuotientAndRemainderResult<Polynomial<T>> qr = polynomials.quotientAndRemainder(t, prevLevel.representative);
		int order = prevLevel.order;
		if (qr.getRemainder().equals(polynomials.zero())) {
			result.add(new Type<>(prevLevel, prevLevel.representative));
			t = polynomials.toUnivariate(qr.getQuotient());
			order--;
		}
		Rationals q = Rationals.q();
		Type<R, RE, RFE>.NewtonPolygon polygon = prevLevel.newtonPolygon(t, order);
		RFE field = prevLevel.reductionFieldExtension.getField();
		for (Fraction slope : polygon.slopes) {
			Fraction lambda = q.negative(slope);
			UnivariatePolynomial<RE> reduced = prevLevel.reduceNextLevel(t, lambda, polygon);
			FactorizationResult<Polynomial<RE>, RE> factorization = field.factorization(reduced);
			for (Polynomial<RE> factor : factorization.primeFactors()) {
				UnivariatePolynomial<RE> nextReduced = field.getUnivariatePolynomialRing().toUnivariate(factor);
				UnivariatePolynomial<T> nextRepresentative = prevLevel.representativeNextLevel(nextReduced, lambda);
				workStack.add(new Type<>(prevLevel, nextRepresentative, nextReduced, lambda,
						factorization.multiplicity(factor), polygon.slopes.size()));
			}
		}
		for (Type<R, RE, RFE> type : workStack) {
			result.addAll(theMontesAlgorithm(t, type));
		}
		return result;
	}

	@Override
	public <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> TheMontesResult<T, S, R, RE, RFE> theMontesAlgorithm(
			UnivariatePolynomial<T> t, Extension<S, R, RE, RFE> reductionAsExtension) {
		if (reductionAsExtension.asVectorMap().evaluate(reductionAsExtension.extension().one()).dimension() != 1) {
			throw new ArithmeticException("Not a trivial extension!");
		}
		List<OkutsuType<T, S, R, RE, RFE>> types = new ArrayList<>();
		UnivariatePolynomial<S> reduced = reduceUnivariatePolynomial(t);
		FactorizationResult<Polynomial<S>, S> reducedFactorization = reduction().factorization(reduced);
		for (Polynomial<S> factor : reducedFactorization.primeFactors()) {
			types.addAll(theMontesAlgorithm(t, new Type<>(t, liftUnivariatePolynomial(factor),
					reducedFactorization.multiplicity(factor), reductionAsExtension)));
		}
		return new TheMontesResult<>(types, t);
	}

	@Override
	public TheMontesResult<T, S, ?, ?, ?> theMontesAlgorithm(UnivariatePolynomial<T> t) {
		return theMontesAlgorithm(t, reduction().getExtension(reduction().getUnivariatePolynomialRing().getVar()));
	}

	@Override
	public <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> UnivariatePolynomial<T> invertInType(
			UnivariatePolynomial<T> t, OkutsuType<T, S, R, RE, RFE> type, int precision) {
		UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
		if (!type.valuation(t).equals(Value.ZERO)) {
			throw new ArithmeticException("Not divisible in ring");
		}
		int requiredDigitAccuracy;
		if (type.quality().isInfinite()) {
			requiredDigitAccuracy = MiscAlgorithms.DivRoundUp(
					2 * type.ramificationIndex()
							* (Math.max(0, precision - type.previousLevel().precision().value()) + 1),
					type.ramificationIndex());
			// ExtendedResultantResult<T> er = polynomials.extendedResultant(t,
			// type.getPolynomial());
//			return er.getCoeff1();
		} else {
			requiredDigitAccuracy = MiscAlgorithms.DivRoundUp(2 * type.ramificationIndex()
					* (Math.max(0, precision - type.previousLevel().precision().value()) + 1) + type.quality().value(),
					type.ramificationIndex()) + 4 * type.exponent();
		}
		int mu = valuationOfUnivariatePolynomial(t).value();
		UnivariatePolynomial<T> denormalized = polynomials.multiply(field.power(field.uniformizer(), -mu), t);
		RE reduced = type.reduction().extension().inverse(type.reduce(denormalized));
		UnivariatePolynomial<T> inverse = polynomials.multiply(field.power(field.uniformizer(), mu),
				type.lift(reduced, -mu * type.ramificationIndex()));
		for (int s = 1; s < precision; s *= 2) {
			inverse = polynomials.toUnivariate(polynomials.remainder(
					polynomials.multiply(inverse,
							polynomials.subtract(polynomials.getInteger(2), polynomials.multiply(t, inverse))),
					type.representative()));
			inverse = roundUnivariatePolynomial(inverse, requiredDigitAccuracy);
		}
		return inverse;
	}

	@Override
	public <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> OkutsuType<T, S, R, RE, RFE> singleFactorLifting(
			OkutsuType<T, S, R, RE, RFE> type, int accuracy) {
		if (type.precision().isInfinite()) {
			return type;
		}
		if (!type.leaf()) {
			throw new ArithmeticException("Not a leaf type!");
		}
		int initialAccuracy = type.previousLevel().precision().value();
		int requiredDigitAccuracy = MiscAlgorithms.DivRoundUp(
				2 * type.ramificationIndex() * (accuracy - initialAccuracy + 1) + type.quality().value(),
				type.ramificationIndex()) + 4 * type.exponent();
		if (field.exactness() != Exactness.EXACT && requiredDigitAccuracy > field.getAccuracy()) {
			return field.withAccuracy(requiredDigitAccuracy).ringOfIntegers().singleFactorLifting(type, accuracy);
		}
		UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
		UnivariatePolynomial<T> f = type.getPolynomial();
		UnivariatePolynomial<T> phi = type.representative();
		QuotientAndRemainderResult<Polynomial<T>> qr = polynomials.quotientAndRemainder(f, phi);
		UnivariatePolynomial<T> a = polynomials.toUnivariate(qr.getQuotient());
		UnivariatePolynomial<T> a0 = polynomials.toUnivariate(qr.getRemainder());
		if (a0.equals(polynomials.zero())) {
			return type;
		}
		UnivariatePolynomial<T> a1 = polynomials.toUnivariate(polynomials.remainder(a, phi));
		int h = type.lambda().getNumerator().getValue().intValueExact();
		UnivariatePolynomial<T> psi = type.lift(type.reduction().extension().one(), -type.valuation(a1).value());
		UnivariatePolynomial<T> capitalA0 = polynomials
				.toUnivariate(polynomials.remainder(polynomials.multiply(a0, psi), phi));
		UnivariatePolynomial<T> capitalA1 = polynomials
				.toUnivariate(polynomials.remainder(polynomials.multiply(a1, psi), phi));
		UnivariatePolynomial<T> a1Inverse = invertInType(capitalA1, type, h);
		UnivariatePolynomial<T> capitalA = polynomials
				.toUnivariate(polynomials.remainder(polynomials.multiply(capitalA0, a1Inverse), phi));
		UnivariatePolynomial<T> capitalPhi = polynomials.add(phi, capitalA);
		UnivariatePolynomial<T> c1Inverse = a1Inverse;
		for (; h < type.ramificationIndex() * (accuracy - initialAccuracy); h *= 2) {
			int digitAccuracy = MiscAlgorithms.DivRoundUp(4 * (h + 1) + type.quality().value(),
					type.ramificationIndex());
			int calculationAccuracy = digitAccuracy + 4 * type.exponent();
			qr = polynomials.quotientAndRemainder(f, capitalPhi);
			UnivariatePolynomial<T> c = polynomials.toUnivariate(qr.getQuotient());
			UnivariatePolynomial<T> c0 = polynomials.toUnivariate(qr.getRemainder());
			UnivariatePolynomial<T> c1 = polynomials.toUnivariate(polynomials.remainder(c, capitalPhi));
			UnivariatePolynomial<T> capitalC0 = polynomials
					.toUnivariate(polynomials.remainder(polynomials.multiply(c0, psi), capitalPhi));
			UnivariatePolynomial<T> capitalC1 = polynomials
					.toUnivariate(polynomials.remainder(polynomials.multiply(c1, psi), capitalPhi));
			c1Inverse = polynomials.toUnivariate(polynomials.remainder(polynomials.multiply(c1Inverse,
					polynomials.subtract(polynomials.getInteger(2), polynomials.multiply(capitalC1, c1Inverse))),
					capitalPhi));
			c1Inverse = roundUnivariatePolynomial(c1Inverse, calculationAccuracy);
			c = polynomials.toUnivariate(polynomials.remainder(polynomials.multiply(c1Inverse, capitalC0), capitalPhi));
			c = roundUnivariatePolynomial(c, digitAccuracy);
			capitalPhi = polynomials.add(capitalPhi, c);
		}
		return type.updateLeaf(capitalPhi);
	}

	// https://arxiv.org/pdf/1507.04058.pdf
	@Override
	public <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> List<List<UnivariatePolynomial<T>>> integralBasis(
			Polynomial<T> minimalPolynomial, TheMontesResult<T, S, R, RE, RFE> theMontes, boolean reduced) {
		if (!theMontes.getPolynomial().equals(minimalPolynomial)) {
			throw new IllegalArgumentException("Not the the Montes result of the given polynomial!");
		}
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
		int numTypes = theMontes.getTypes().size();
		// TODO replace with explicit formulas
		List<List<UnivariatePolynomial<T>>> result = new ArrayList<>();
		List<Set<Integer>> improvements = new ArrayList<>();
		if (numTypes > 1) {
			List<List<Fraction>> phiValuations = new ArrayList<>();
			List<Fraction> differences = new ArrayList<>();
			for (int i = 0; i < numTypes; i++) {
				OkutsuType<T, S, R, RE, RFE> type = theMontes.getTypes().get(i);
				differences.add(q.zero());
				List<Fraction> phiValuation = new ArrayList<>();
				phiValuations.add(phiValuation);
				for (int k = 0; k < numTypes; k++) {
					phiValuation.add(scaledValue(type, theMontes.getTypes().get(k).representative()));
				}
			}
			for (int i = 0; i < numTypes; i++) {
				OkutsuType<T, S, R, RE, RFE> type = theMontes.getTypes().get(i);
				Set<Integer> improved = new TreeSet<>();
				improvements.add(improved);
				boolean improvedStable = false;
				while (!improvedStable) {
					improvedStable = true;
					for (int k = 0; k < numTypes; k++) {
						if (i == k) {
							continue;
						}
						OkutsuType<T, S, R, RE, RFE> secondaryType = theMontes.getTypes().get(k);
						Fraction difference = q.zero();
						boolean requiresImprovement = false;
						for (int j = 0; j < type.representative().degree(); j++) {
							UnivariatePolynomial<T> integralGenerator = type.divisorPolynomial(j);
							Fraction generatorValue = scaledValue(type, integralGenerator);
							for (int improvedIndex : improved) {
								generatorValue = q.add(phiValuations.get(i).get(improvedIndex), generatorValue);
							}
							Fraction secondaryValue = scaledValue(secondaryType, integralGenerator);
							for (int improvedIndex : improved) {
								secondaryValue = q.add(phiValuations.get(k).get(improvedIndex), secondaryValue);
							}
							requiresImprovement = requiresImprovement
									|| requiresImprovement(i, k, generatorValue, secondaryValue, reduced);
							difference = q.maximum(
									q.subtract(generatorValue, scaledValue(secondaryType, integralGenerator)),
									difference);
						}
						if (requiresImprovement && !improved.contains(k)) {
							improved.add(k);
							improvedStable = false;
						} else if (requiresImprovement) {
							differences.set(k, q.maximum(differences.get(k), difference));
						}
					}
				}
			}
			for (int i = 0; i < numTypes; i++) {
				OkutsuType<T, S, R, RE, RFE> type = theMontes.getTypes().get(i);
				Fraction difference = differences.get(i);
				if (!reduced) {
					difference = q.getInteger(difference.roundUp());
				}
				if (difference.equals(q.zero())) {
					continue;
				}
				singleFactorLifting(type,
						q.multiply(difference, q.getInteger(type.ramificationIndex())).asInteger().intValueExact());
			}
		} else {
			improvements.add(Collections.emptySet());
		}
		for (int i = 0; i < numTypes; i++) {
			OkutsuType<T, S, R, RE, RFE> type = theMontes.getTypes().get(i);
			List<UnivariatePolynomial<T>> typeResult = new ArrayList<>();
			result.add(typeResult);
			UnivariatePolynomial<T> multiplier = polynomials.one();
			for (int improve : improvements.get(i)) {
				multiplier = polynomials.multiply(theMontes.getTypes().get(improve).representative(), multiplier);
			}
			for (int j = 0; j < type.representative().degree(); j++) {
				UnivariatePolynomial<T> integralGenerator = type.divisorPolynomial(j);
				integralGenerator = polynomials.multiply(multiplier, integralGenerator);
				Fraction generatorValue = scaledValue(type, integralGenerator);
				typeResult.add(roundUnivariatePolynomial(
						polynomials.multiply(field.power(field.uniformizer(),
								-MiscAlgorithms.DivRoundDown(generatorValue.numeratorIntValueExact(),
										generatorValue.denominatorIntValueExact())),
								integralGenerator),
						reduced ? 2 : 1));
			}
		}
		return result;
	}

	private <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> Fraction scaledValue(
			OkutsuType<T, S, R, RE, RFE> type, UnivariatePolynomial<T> t) {
		Rationals q = Rationals.q();
		return q.getFraction(type.valuation(t).value(), type.ramificationIndex());
	}

	private boolean requiresImprovement(int i, int k, Fraction firstValue, Fraction secondValue, boolean reduced) {
		if (reduced) {
			return firstValue.compareTo(secondValue) >= 0;
		}
		int cmp = firstValue.roundDown().compareTo(secondValue.roundDown());
		if (i < k) {
			return cmp > 0;
		}
		return cmp >= 0;
	}

	@Override
	public List<UnivariatePolynomial<T>> triagonalizeIntegralBasis(Polynomial<T> minimalPolynomial,
			List<List<UnivariatePolynomial<T>>> integralBasis) {
		List<UnivariatePolynomial<T>> flattenedBasis = new ArrayList<>();
		for (List<UnivariatePolynomial<T>> typeResult : integralBasis) {
			flattenedBasis.addAll(typeResult);
		}
//		UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
//		List<UnivariatePolynomial<T>> basis = new ArrayList<>();
//		Collections.sort(flattenedBasis, new Comparator<>() {
//			@Override
//			public int compare(UnivariatePolynomial<T> o1, UnivariatePolynomial<T> o2) {
//				if (o1.degree() != o2.degree()) {
//					return o1.degree() - o2.degree();
//				}
//				return field.valuation(o1.leadingCoefficient()).compareTo(field.valuation(o2.leadingCoefficient()));
//			}
//		});
//		for (int i = 0; i < minimalPolynomial.degree(); i++) {
//			UnivariatePolynomial<T> polynomial = flattenedBasis.get(i);
//			T coeff = polynomial.univariateCoefficient(i);
//			Value pivotValue = field.valuation(coeff);
//			T multiplier = field.multiply(field.inverse(coeff), field.power(uniformizer(), pivotValue.value()));
//			basis.add(roundUnivariatePolynomial(polynomials.scalarMultiply(multiplier, polynomial), 1));
//		}
//		return basis;
		UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
		List<UnivariatePolynomial<T>> basis = new ArrayList<>();
		for (int i = minimalPolynomial.degree() - 1; i >= 0; i--) {
			Value minValue = Value.INFINITY;
			int minValueIndex = -1;
			for (int j = 0; j < flattenedBasis.size(); j++) {
				Value value = field.valuation(flattenedBasis.get(j).univariateCoefficient(i));
				if (value.compareTo(minValue) < 0) {
					minValue = value;
					minValueIndex = j;
				}
			}
			UnivariatePolynomial<T> pivot = flattenedBasis.remove(minValueIndex);
			T pivotCoeff = pivot.univariateCoefficient(i);
			for (int j = 0; j < flattenedBasis.size(); j++) {
				UnivariatePolynomial<T> generator = flattenedBasis.get(j);
				T generatorCoeff = generator.univariateCoefficient(i);
				if (generatorCoeff.equals(zero())) {
					continue;
				}
				T shear = divideChecked(negative(generatorCoeff), pivotCoeff);
				flattenedBasis.set(j, polynomials.add(generator, polynomials.scalarMultiply(shear, pivot)));
			}
			Value pivotValue = field.valuation(pivotCoeff);
			T multiplier = field.multiply(field.inverse(pivotCoeff), field.power(uniformizer(), pivotValue.value()));
			basis.add(roundUnivariatePolynomial(polynomials.scalarMultiply(multiplier, pivot), 1));
		}
		Collections.reverse(basis);
		return basis;

//		FreeModule<T> vectors = new FreeModule<>(this, minimalPolynomial.degree());
//		Value minValue = Value.INFINITY;
//		for (int i = 0; i < minimalPolynomial.degree(); i++) {
//			minValue = minValue.min(valuationOfUnivariatePolynomial(flattenedBasis.get(i)));
//		}
//		T power = power(uniformizer(), -minValue.value());
//		T inversePower = field.inverse(power);
//		List<Vector<T>> asVectors = new ArrayList<>();
//		for (int i = 0; i < minimalPolynomial.degree(); i++) {
//			Vector<T> reverse = polynomials.asVector(
//					polynomials.scalarMultiply(power, flattenedBasis.get(minimalPolynomial.degree() - i - 1)),
//					minimalPolynomial.degree() - 1);
//			List<T> asList = new ArrayList<>();
//			asList.addAll(reverse.asList());
//			Collections.reverse(asList);
//			asVectors.add(new Vector<>(asList));
//		}
//		Matrix<T> m = Matrix.fromRows(asVectors);
//		// MatrixModule<T>.LDUPResult lup = vectors.matrixAlgebra().ldup(m);
//		// if (!lup.getDeterminant().equals(field.one())) {
//		// throw new ArithmeticException("Triangulation did not work as expected");
//		// }
//		Matrix<T> triangBasis = vectors.matrixAlgebra().lupUpperTriangle(m);
//		List<UnivariatePolynomial<T>> basis = new ArrayList<>();
//		for (int i = 0; i < minimalPolynomial.degree(); i++) {
//			Vector<T> reverse = triangBasis.row(minimalPolynomial.degree() - i);
//			List<T> asList = new ArrayList<>();
//			asList.addAll(reverse.asList());
//			Collections.reverse(asList);
//			UnivariatePolynomial<T> basisPolynomial = polynomials.multiply(inversePower,
//					polynomials.getPolynomial(asList));
//			T lc = basisPolynomial.leadingCoefficient();
//			int lcValue = valuation(lc).value();
//			lc = field.multiply(lc, field.power(field.uniformizer(), -lcValue));
//			T inv = inverse(lc);
//			basisPolynomial = roundUnivariatePolynomial(polynomials.multiply(inv, basisPolynomial), 1);
//			basis.add(basisPolynomial);
//		}
//		return basis;
	}

	@Override
	public FactorizationResult<Polynomial<T>, T> factorization(UnivariatePolynomial<T> t) {
		if (isComplete()) {
			return factorization(t, false, field.getAccuracy());
		}
		UnivariatePolynomialRing<T> polynomials = getUnivariatePolynomialRing();
		T content = polynomials.content(t);
		FactorizationResult<Polynomial<T>, T> factorization = field.factorization(polynomials.contentFree(t));
		SortedMap<Polynomial<T>, Integer> result = new TreeMap<>();
		for (Polynomial<T> factor : factorization.primeFactors()) {
			Value val = Value.INFINITY;
			UnivariatePolynomial<T> p = polynomials.toUnivariate(factor);
			for (int i = 0; i <= p.degree(); i++) {
				val = val.min(valuation(p.univariateCoefficient(i)));
			}
			if (val.compareTo(Value.ZERO) < 0) {
				T uniformizerPower = power(uniformizer(), -val.value());
				p = polynomials.multiply(uniformizerPower, p);
				content = divide(content, uniformizerPower);
			}
			result.put(p, factorization.multiplicity(factor));
		}
		Value val = valuation(content);
		if (val.value() != 0) {
			result.put(polynomials.getEmbedding(uniformizer()), val.value());
			content = divide(content, power(uniformizer(), val.value()));
		}
		return new FactorizationResult<>(content, result);
	}

	@Override
	public FactorizationResult<Polynomial<T>, T> factorization(UnivariatePolynomial<T> t, int accuracy) {
		return factorization(t, false, accuracy);
	}

	@Override
	public FactorizationResult<Polynomial<T>, T> factorization(UnivariatePolynomial<T> t, boolean forField,
			int accuracy) {
		if (!isComplete()) {
			throw new ArithmeticException("Factorization only works for complete local rings!");
		}
		UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
		if (t.equals(polynomials.zero())) {
			throw new ArithmeticException("Cannot factor zero!");
		}
		SortedMap<Polynomial<T>, Integer> result = new TreeMap<>();
		T unit;
		if (forField) {
			unit = t.leadingCoefficient();
			Value lcValue = field.valuation(t.leadingCoefficient());
			Value minValue = Value.ZERO;
			ValueGroup g = ValueGroup.g();
			for (int i = t.degree(); i >= 0; i--) {
				Value v = g.subtract(field.valuation(t.univariateCoefficient(i)), lcValue);
				v = g.multiply(t.degree() - i, v);
				minValue = g.min(minValue, v);
			}
			int power = 0;
			if (minValue.compareTo(Value.ZERO) < 0) {
				int value = -minValue.value();
				power = MiscAlgorithms.DivRoundUp(value, t.degree());
			}
			t = polynomials.normalize(polynomials.substitute(t,
					Collections.singletonList(polynomials.getEmbedding(field.power(field.uniformizer(), -power), 1))));
			FactorizationResult<Polynomial<T>, T> factors = factorization(t, field.exact(), accuracy);
			unit = multiply(factors.getUnit(), unit);
			for (Polynomial<T> factor : factors.primeFactors()) {
				result.put(
						polynomials.normalize(polynomials.substitute(factor,
								Collections.singletonList(
										polynomials.getEmbedding(field.power(field.uniformizer(), power), 1)))),
						factors.multiplicity(factor));
			}

		} else {
			FactorizationResult<T, T> content = uniqueFactorization(getUnivariatePolynomialRing().content(t));
			unit = content.getUnit();
			for (T contentFactor : content.primeFactors()) {
				result.put(polynomials.getEmbedding(contentFactor), content.multiplicity(contentFactor));
			}
			FactorizationResult<Polynomial<T>, T> factors = factorization(t, field.exact(), accuracy);
			unit = multiply(factors.getUnit(), unit);
			result.putAll(factors.factorMap());
		}
		return new FactorizationResult<>(unit, result);
	}

	private <E extends Element<E>, LE extends DiscreteValuationField<E, S>> FactorizationResult<Polynomial<T>, T> factorization(
			UnivariatePolynomial<T> t, OtherVersion<T, E, S, LE> exact, int accuracy) {
		LE exactField = exact.getField();
		DiscreteValuationRing<E, S> exactRing = exactField.ringOfIntegers();
		UnivariatePolynomialRing<E> exactPolynomials = exactField.getUnivariatePolynomialRing();
		UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
		FactorizationResult<Polynomial<E>, E> exactSquareFreeFactors = exactPolynomials
				.squareFreeFactorization(exactPolynomials.getEmbedding(t, exact.getRetraction()));
		SortedMap<Polynomial<T>, Integer> result = new TreeMap<>();
		for (Polynomial<E> squareFree : exactSquareFreeFactors.primeFactors()) {
			List<Polynomial<T>> factors = factorSquareFreePolynomial(squareFree,
					theMontesAlgorithm(polynomials.getEmbedding(squareFree, exact.getEmbedding())), exactRing,
					accuracy);
			for (Polynomial<T> factor : factors) {
				result.put(roundPolynomial(factor, accuracy), exactSquareFreeFactors.multiplicity(squareFree));
			}
		}
		return new FactorizationResult<>(exact.getEmbedding().evaluate(exactSquareFreeFactors.getUnit()), result);
	}

	private <E extends Element<E>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> List<Polynomial<T>> factorSquareFreePolynomial(
			Polynomial<E> squareFree, TheMontesResult<T, S, R, RE, RFE> theMontes,
			DiscreteValuationRing<E, S> exactRing, int accuracy) {
		List<Polynomial<T>> result = new ArrayList<>();
		for (OkutsuType<T, S, R, RE, RFE> type : theMontes.getTypes()) {
			result.add(singleFactorLifting(type, accuracy).representative());
		}
		return result;
	}

	@Override
	public LocalFieldAsModule<T, S> fieldOfFractionsAsModule() {
		return new LocalFieldAsModule<>(field);
	}

	public static class LocalFieldAsModule<T extends Element<T>, S extends Element<S>> extends AbstractModule<T, T>
			implements Module<T, T> {
		private DiscreteValuationField<T, S> field;

		private LocalFieldAsModule(DiscreteValuationField<T, S> field) {
			this.field = field;
		}

		public DiscreteValuationField<T, S> getField() {
			return field;
		}

		@Override
		public DiscreteValuationRing<T, S> getRing() {
			return field.ringOfIntegers();
		}

		@Override
		public T scalarMultiply(T t, T s) {
			return field.multiply(t, s);
		}

		@Override
		public boolean isFree() {
			return true;
		}

		@Override
		public Ideal<T> annihilator() {
			return getRing().getZeroIdeal();
		}
		
		@Override
		public List<Vector<T>> getSyzygies() {
			return Collections.emptyList();
		}

		@Override
		public boolean isLinearIndependent(List<T> s) {
			List<T> nonzero = new ArrayList<>();
			for (T t : s) {
				if (!t.equals(zero())) {
					nonzero.add(t);
				}
			}
			if (nonzero.size() < 2) {
				return true;
			}
			return false;
		}

		@Override
		public boolean isGeneratingModule(List<T> s) {
			return false;
		}

		@Override
		public List<Vector<T>> nonTrivialCombinations(List<T> s) {
			int firstValue = -1;
			int firstIndex = -1;
			List<Vector<T>> combinations = new ArrayList<>();
			for (int i = 0; i < s.size(); i++) {
				Value v = field.valuation(s.get(i));
				if (v.isInfinite()) {
					continue;
				}
				int value = v.value();
				if (firstIndex == -1) {
					firstValue = value;
					firstIndex = i;
					continue;
				}
				int minValue = Math.min(firstValue, value);
				minValue = Math.min(minValue, 0);
				T power = field.power(field.uniformizer(), -minValue);
				T first = field.multiply(power, s.get(firstIndex));
				T current = field.multiply(power, s.get(i));
				ExtendedEuclideanResult<T> ee = field.ringOfIntegers().extendedEuclidean(first, current);
				List<T> combination = new ArrayList<>();
				for (int j = 0; j < firstIndex; j++) {
					combination.add(zero());
				}
				combination.add(ee.getCoeff1());
				for (int j = firstIndex + 1; j < i; j++) {
					combination.add(zero());
				}
				combination.add(ee.getCoeff2());
				for (int j = i + 1; j < s.size(); j++) {
					combination.add(zero());
				}
				combinations.add(new Vector<>(combination));
			}
			return combinations;
		}

		@Override
		public List<T> getModuleGenerators() {
			throw new InfinityException();
		}

		@Override
		public Vector<T> asVector(T s) {
			Value value = field.valuation(s);
			if (value.compareTo(Value.ZERO) >= 0) {
				return new Vector<>(Collections.singletonList(s));
			}
			int power = -value.value();
			List<T> result = new ArrayList<>();
			for (int i = 0; i < power; i++) {
				result.add(zero());
			}
			result.add(field.multiply(field.power(field.uniformizer(), power), s));
			return new Vector<>(result);
		}

		@Override
		public T fromVector(Vector<T> vector) {
			T result = zero();
			for (int i = 0; i < vector.dimension(); i++) {
				result = add(result, field.multiply(field.power(field.uniformizer(), i), vector.get(i + 1)));
			}
			return result;
		}

		@Override
		public T zero() {
			return field.zero();
		}

		@Override
		public T add(T s1, T s2) {
			return field.add(s1, s2);
		}

		@Override
		public T negative(T s) {
			return field.negative(s);
		}

		@Override
		public Exactness exactness() {
			return field.exactness();
		}

		@Override
		public T getRandomElement() {
			return field.getRandomElement();
		}

		@Override
		public boolean isFinite() {
			return field.isFinite();
		}

		@Override
		public BigInteger getNumberOfElements() throws InfinityException {
			return field.getNumberOfElements();
		}

		@Override
		public Iterator<T> iterator() {
			return field.iterator();
		}

	}

	private static class LocalIdeal<T extends Element<T>, S extends Element<S>> extends AbstractIdeal<T> {
		private Value valuation;
		private T generator;
		private DiscreteValuationRing<T, S> ring;

		private LocalIdeal(DiscreteValuationRing<T, S> ring, Value valuation) {
			super(ring);
			this.ring = ring;
			this.valuation = valuation;
			this.generator = valuation.isInfinite() ? ring.zero() : ring.power(ring.uniformizer(), valuation.value());
		}

		@Override
		public boolean isPrimary() {
			return !valuation.isInfinite();
		}

		@Override
		public boolean isPrime() {
			return valuation.isInfinite() || valuation.value() == 1;
		}

		@Override
		public boolean isMaximal() {
			return valuation.equals(new Value(1));
		}

//		@Override
//		public List<List<T>> nonTrivialCombinations(List<T> s) {
//			Matrix<T> m = new Matrix<>(Collections.singletonList(s));
//			MatrixModule<T> mm = new MatrixModule<>(ring, 1, s.size());
//			List<Vector<T>> kernelBasis = mm.kernelBasis(m);
//			List<List<T>> result = new ArrayList<>();
//			for (Vector<T> basisVector : kernelBasis) {
//				result.add(basisVector.asList());
//			}
//			return result;
//		}

		@Override
		public List<T> generators() {
			return Collections.singletonList(generator);
		}

		@Override
		public List<T> generate(T t) {
			if (valuation.isInfinite()) {
				return Collections.singletonList(ring.zero());
			}
			return Collections.singletonList(ring.divide(ring.subtract(t, residue(t)), generator));
		}

		@Override
		public T residue(T t) {
			if (valuation.isInfinite()) {
				return t;
			}
			return ring.round(t, valuation.value());
//			T result = ring.zero();
//			for (int i = 0; i < valuation; i++) {
//				T level = ring.power(ring.uniformizer(), i);
//				if (level.equals(zero())) {
//					break;
//				}
//				result = ring.add(result,
//						ring.multiply(ring.lift(ring.reduce(ring.divide(ring.subtract(t, result), level))), level));
//			}
//			return result;
		}

		@Override
		public boolean contains(T t) {
			return ring.valuation(t).compareTo(valuation) >= 0;
		}

		@Override
		public Value maximumPowerContains(T t) {
			if (t.equals(zero()) || valuation.equals(Value.ZERO)) {
				return Value.INFINITY;
			}
			if (valuation.isInfinite()) {
				return Value.ZERO;
			}
			return new Value(ring.valuation(t).value() / valuation.value());
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
		public String toString() {
			return "(" + generator + ")";
		}
	}
}
