package fields.polynomials;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.helper.AbstractElement;
import fields.helper.AbstractIdeal;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.local.FormalPowerSeries;
import fields.local.FormalPowerSeries.PowerSeries;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.LocalizedCoordinateRing.LocalizedElement;
import fields.vectors.Vector;
import varieties.FunctionField;
import varieties.RationalFunction;
import varieties.projective.GenericProjectiveScheme;

public class CoordinateRing<T extends Element<T>> extends AbstractAlgebra<Polynomial<T>, CoordinateRingElement<T>>
		implements Ring<CoordinateRingElement<T>> {
	private PolynomialRing<T> ring;
	private PolynomialIdeal<T> ideal;
	private int dimension;
	private int degree;
	private Set<Integer> boundVariables;

	public static class CoordinateRingElement<T extends Element<T>> extends AbstractElement<CoordinateRingElement<T>> {
		private Polynomial<T> polynomial;

		public CoordinateRingElement(Ideal<Polynomial<T>> ideal, Polynomial<T> polynomial) {
			this.polynomial = ideal.residue(polynomial);
		}

		private CoordinateRingElement(Ideal<Polynomial<T>> ideal, Polynomial<T> polynomial, boolean reduce) {
			if (reduce) {
				this.polynomial = ideal.residue(polynomial);
			} else {
				this.polynomial = polynomial;
			}
		}

		public String toString() {
			return this.polynomial.toString();
		}

		public Polynomial<T> getElement() {
			return this.polynomial;
		}

		@Override
		public int compareTo(CoordinateRingElement<T> o) {
			return this.polynomial.compareTo(o.polynomial);
		}
	}

	public static <T extends Element<T>> CoordinateRing<T> tensorProduct(CoordinateRing<T> t1, CoordinateRing<T> t2) {
		return tensorProduct(t1, t2, Monomial.GREVLEX);
	}

	public static <T extends Element<T>> CoordinateRing<T> tensorProduct(CoordinateRing<T> t1, CoordinateRing<T> t2,
			Comparator<Monomial> comparator) {
		if (!t1.getPolynomialRing().getRing().equals(t2.getPolynomialRing().getRing())) {
			throw new ArithmeticException("Not over common base!");
		}
		Ring<T> base = t1.getPolynomialRing().getRing();
		PolynomialRing<T> polynomialRing = AbstractPolynomialRing.getPolynomialRing(base,
				t1.getPolynomialRing().numberOfVariables() + t2.getPolynomialRing().numberOfVariables(), comparator);
		int[] map1 = new int[t1.getPolynomialRing().numberOfVariables()];
		for (int i = 0; i < map1.length; i++) {
			map1[i] = i;
		}
		int[] map2 = new int[t2.getPolynomialRing().numberOfVariables()];
		for (int i = 0; i < map2.length; i++) {
			map2[i] = i + map1.length;
		}
		List<Polynomial<T>> polynomials = new ArrayList<>();
		for (Polynomial<T> polynomial : t1.getIdeal().generators()) {
			polynomials.add(polynomialRing.getEmbedding(polynomial, map1));
		}
		for (Polynomial<T> polynomial : t2.getIdeal().generators()) {
			polynomials.add(polynomialRing.getEmbedding(polynomial, map2));
		}
		return polynomialRing.getIdeal(polynomials).divideOut();
	}

	public static <T extends Element<T>> CoordinateRing<T> getCoordinateRing(CoordinateIdeal<T> ideal) {
		return getCoordinateRing(ideal.asPolynomialIdeal());
	}

	public static <T extends Element<T>> CoordinateRing<T> getCoordinateRing(PolynomialIdeal<T> ideal) {
		return ideal.divideOut();
	}

	CoordinateRing(PolynomialRing<T> ring, PolynomialIdeal<T> ideal) {
		this.ring = ring;
		this.ideal = ideal;
		int i = 0;
		int nonTrivial = 0;
		for (Polynomial<T> generator : this.ideal.generators()) {
			if (generator.degree() > 0) {
				nonTrivial++;
			}
		}
		int[][] leadingMonomials = new int[nonTrivial][this.ring.numberOfVariables()];
		for (Polynomial<T> generator : this.ideal.generators()) {
			if (generator.degree() <= 0) {
				continue;
			}
			leadingMonomials[i] = generator.leadingMonomial().exponents();
			i++;
		}
		this.boundVariables = new TreeSet<>();
		this.degree = -1;
		if (ideal.contains(ring.one())) {
			this.dimension = -1;
		} else {
			int[] set = this.makeset(0, new int[this.ring.numberOfVariables() + 1], leadingMonomials);
			this.dimension = this.ring.krullDimension() - set[this.ring.numberOfVariables()];
			for (int j = 0; j < ring.numberOfVariables(); j++) {
				if (set[j] != 0) {
					boundVariables.add(j + 1);
				}
			}
		}
	}

	private int[] makeset(int level, int[] set, int[][] leadingMonomials) {
		if (level == leadingMonomials.length)
			return set;
		int[] optimalset = new int[set.length];
		optimalset[set.length - 1] = -1;
		for (int i = 0; i < leadingMonomials[level].length; i++) {
			if (leadingMonomials[level][i] == 0)
				continue;
			if (set[i] == 1) {
				return this.makeset(level + 1, set, leadingMonomials);
			}
			set[i] = 1;
			set[set.length - 1]++;
			int[] sethere = this.makeset(level + 1, set, leadingMonomials);
			if (optimalset[set.length - 1] == -1 || optimalset[set.length - 1] > sethere[set.length - 1]) {
				optimalset = Arrays.copyOf(sethere, sethere.length);
			}
			set[i] = 0;
			set[set.length - 1]--;
		}
		return optimalset;
	}

	@Override
	public Exactness exactness() {
		return ring.exactness();
	}

	@Override
	public String toString() {
		return this.ring.toString() + "/" + this.ideal.toString();
	}

	public PolynomialRing<T> getPolynomialRing() {
		return ring;
	}

	public Ring<Polynomial<T>> getRing() {
		return ring;
	}

	public PolynomialIdeal<T> getIdeal() {
		return ideal;
	}

	@Override
	public int krullDimension() {
		return dimension;
	}

	@Override
	public List<Ideal<CoordinateRingElement<T>>> maximalPrimeIdealChain(Ideal<CoordinateRingElement<T>> start) {
		List<Ideal<Polynomial<T>>> polynomialChain = ring
				.maximalPrimeIdealChain(((CoordinateIdeal<T>) start).asPolynomialIdeal());
		List<Ideal<CoordinateRingElement<T>>> result = new ArrayList<>();
		for (Ideal<Polynomial<T>> polynomialIdeal : polynomialChain) {
			result.add(getIdeal(polynomialIdeal));
		}
		return result;
	}

	@Override
	public List<Ideal<CoordinateRingElement<T>>> maximalPrimeIdealChain(Ideal<CoordinateRingElement<T>> start,
			Ideal<CoordinateRingElement<T>> end) {
		List<Ideal<Polynomial<T>>> polynomialChain = ring.maximalPrimeIdealChain(
				((CoordinateIdeal<T>) start).asPolynomialIdeal(), ((CoordinateIdeal<T>) end).asPolynomialIdeal());
		List<Ideal<CoordinateRingElement<T>>> result = new ArrayList<>();
		for (Ideal<Polynomial<T>> polynomialIdeal : polynomialChain) {
			result.add(getIdeal(polynomialIdeal));
		}
		return result;
	}

	public Set<Integer> boundVariables() {
		return boundVariables;
	}

	private <S extends Element<S>> FieldOfFractionsResult<CoordinateRingElement<T>, RationalFunction<S>> fieldOfFractions(
			FieldOfFractionsResult<T, S> baseFractions) {
		PolynomialRing<S> polynomialRing = AbstractPolynomialRing.getPolynomialRing(baseFractions.getField(),
				ring.numberOfVariables(), ring.getComparator());
		PolynomialIdeal<S> fractionIdeal = polynomialRing.getEmbedding(ideal, baseFractions.getEmbedding());
		FunctionField<S> functionField = new FunctionField<>(
				GenericProjectiveScheme.fromAffineCoordinateRing(fractionIdeal.divideOut()));
		MathMap<RationalFunction<S>, S> denominator = new MathMap<>() {

			@Override
			public S evaluate(RationalFunction<S> t) {
				Polynomial<S> numerator = t.getNumerator();
				T denominator = baseFractions.getRing().one();
				for (Monomial m : numerator.monomials()) {
					denominator = baseFractions.getRing().lcm(denominator,
							baseFractions.getDenominator().evaluate(numerator.coefficient(m)));
				}
				Polynomial<S> denom = t.getDenominator();
				for (Monomial m : denom.monomials()) {
					denominator = baseFractions.getRing().lcm(denominator,
							baseFractions.getDenominator().evaluate(denom.coefficient(m)));
				}
				return baseFractions.getEmbedding().evaluate(denominator);
			}
		};
		return new FieldOfFractionsResult<>(this, functionField, new MathMap<>() {

			@Override
			public RationalFunction<S> evaluate(CoordinateRingElement<T> t) {
				return functionField
						.getEmbedding(polynomialRing.getEmbedding(t.getElement(), baseFractions.getEmbedding()));
			}
		}, new MathMap<>() {

			@Override
			public CoordinateRingElement<T> evaluate(RationalFunction<S> t) {
				return getEmbedding(
						ring.getEmbedding(polynomialRing.multiply(denominator.evaluate(t), t.getNumerator()),
								baseFractions.getAsInteger()));
			}
		}, new MathMap<>() {

			@Override
			public CoordinateRingElement<T> evaluate(RationalFunction<S> t) {
				return getEmbedding(
						ring.getEmbedding(polynomialRing.multiply(denominator.evaluate(t), t.getDenominator()),
								baseFractions.getAsInteger()));
			}
		}, new MathMap<>() {

			@Override
			public CoordinateRingElement<T> evaluate(RationalFunction<S> t) {
				return getEmbedding(
						ring.getEmbedding(polynomialRing.divideChecked(t.getNumerator(), t.getDenominator()),
								baseFractions.getAsInteger()));
			}
		});
	}

	@Override
	public FieldOfFractionsResult<CoordinateRingElement<T>, ?> fieldOfFractions() {
		return fieldOfFractions(getPolynomialRing().getRing().fieldOfFractions());
	}

	@Override
	public LocalizeResult<CoordinateRingElement<T>, ?, ?, ?> localizeAtIdeal(
			Ideal<CoordinateRingElement<T>> primeIdeal) {
		CoordinateIdeal<T> ideal = (CoordinateIdeal<T>) primeIdeal;
		return localizeAtIdeal(ideal, ring.getRing().localizeAtIdeal(ideal.asPolynomialIdeal().intersectToRing()));
	}

//TODO: FIXME
	private <S extends Element<S>, U extends Element<U>, R extends Element<R>> LocalizeResult<CoordinateRingElement<T>, LocalizedElement<S>, LocalizedElement<S>, S> localizeAtIdeal(
			CoordinateIdeal<T> primeIdeal, LocalizeResult<T, S, U, R> localizeRing) {
		PolynomialRing<S> polynomialRing = AbstractPolynomialRing.getPolynomialRing(localizeRing.getLocalizedRing(),
				ring.numberOfVariables(), ring.getComparator());
		List<Polynomial<S>> mainGenerators = new ArrayList<>();
		for (Polynomial<T> generator : ideal.generators()) {
			mainGenerators.add(polynomialRing.getEmbedding(generator, localizeRing.getEmbedding()));
		}
		PolynomialIdeal<S> localizedMainIdeal = polynomialRing.getIdeal(mainGenerators);
		CoordinateRing<S> coordinateRing = localizedMainIdeal.divideOut();
		List<Polynomial<S>> generators = new ArrayList<>();
		for (CoordinateRingElement<T> generator : primeIdeal.generators()) {
			if (generator.getElement().degree() > 0) {
				generators.add(polynomialRing.getEmbedding(generator.getElement(), localizeRing.getEmbedding()));
			}
		}
		PolynomialIdeal<S> localizedIdeal = polynomialRing.getIdeal(generators);
		LocalizedCoordinateRing<S> result = new LocalizedCoordinateRing<>((Field<S>) localizeRing.getLocalizedRing(),
				coordinateRing, coordinateRing.getIdeal(localizedIdeal));
		MathMap<LocalizedElement<S>, S> denominator = new MathMap<>() {

			@Override
			public S evaluate(LocalizedElement<S> t) {
				Polynomial<S> numerator = t.asPolynomialFraction().getNumerator();
				T denominator = ring.getRing().one();
				for (Monomial m : numerator.monomials()) {
					denominator = ring.getRing().lcm(denominator,
							localizeRing.getDenominator().evaluate(numerator.coefficient(m)));
				}
				Polynomial<S> denom = t.asPolynomialFraction().getDenominator();
				for (Monomial m : denom.monomials()) {
					denominator = ring.getRing().lcm(denominator,
							localizeRing.getDenominator().evaluate(denom.coefficient(m)));
				}
				return localizeRing.getEmbedding().evaluate(denominator);
			}
		};
		return new LocalizeResult<>(this, primeIdeal, result.ringOfIntegers(), new MathMap<>() {

			@Override
			public LocalizedElement<S> evaluate(CoordinateRingElement<T> t) {
				return result.getEmbedding(polynomialRing.getEmbedding(t.getElement(), localizeRing.getEmbedding()));
			}
		}, new MathMap<>() {

			@Override
			public CoordinateRingElement<T> evaluate(LocalizedElement<S> t) {
				Polynomial<S> numerator = t.asPolynomialFraction().getNumerator();
				numerator = polynomialRing.multiply(denominator.evaluate(t), numerator);
				return getEmbedding(ring.getEmbedding(numerator, localizeRing.getAsInteger()));
			}
		}, new MathMap<>() {

			@Override
			public CoordinateRingElement<T> evaluate(LocalizedElement<S> t) {
				Polynomial<S> denom = t.asPolynomialFraction().getDenominator();
				denom = polynomialRing.multiply(denominator.evaluate(t), denom);
				return getEmbedding(ring.getEmbedding(denom, localizeRing.getAsInteger()));
			}
		}, new MathMap<>() {

			@Override
			public CoordinateRingElement<T> evaluate(LocalizedElement<S> t) {
				return getEmbedding(
						ring.getEmbedding(t.asPolynomialFraction().asInteger(), localizeRing.getAsInteger()));
			}
		});
	}

	private static IntE hilbertFunction(int degree, int variable, int numberOfVariables, Set<Monomial> monomials) {
		Integers z = Integers.z();
		if (variable > numberOfVariables) {
			return z.zero();
		}
		if (monomials.size() == 0) {
			return Integers.z().binomialCoefficient(degree + numberOfVariables - variable,
					numberOfVariables - variable);
		}
		Set<Monomial> xMinusOneMonomials = new TreeSet<>();
		Set<Monomial> yMonomials = new TreeSet<>();
		for (Monomial m : monomials) {
			if (m.degree() == 0) {
				return z.zero();
			}
			if (m.exponents()[variable - 1] > 0) {
				int[] exponents = Arrays.copyOf(m.exponents(), numberOfVariables);
				exponents[variable - 1]--;
				xMinusOneMonomials.add(new Monomial(Monomial.GREVLEX, exponents));
			} else {
				if (m.degree() < degree) {
					xMinusOneMonomials.add(m);
				}
				yMonomials.add(m);
			}
		}
		return z.add(hilbertFunction(degree - 1, variable, numberOfVariables, xMinusOneMonomials),
				hilbertFunction(degree, variable + 1, numberOfVariables, yMonomials));
	}

	public IntE hilbertFunction(int degree) {
		Set<Monomial> monomials = new TreeSet<>();
		for (Polynomial<T> generator : ideal.generators()) {
			generator = ring.homogenize(generator);
			if (generator.leadingMonomial().degree() <= degree) {
				monomials.add(generator.leadingMonomial());
			}
		}
		return hilbertFunction(degree, 1, ring.numberOfVariables() + 1, monomials);
	}

	public PowerSeries<Fraction> hilbertSeries(FormalPowerSeries<Fraction> powerSeries) {
		Rationals q = Rationals.q();
		List<Fraction> coefficients = new ArrayList<>();
		for (int i = 0; i <= powerSeries.getAccuracy(); i++) {
			coefficients.add(q.getInteger(hilbertFunction(i)));
		}
		UnivariatePolynomial<Fraction> result = q.getUnivariatePolynomialRing().getPolynomial(coefficients);
		return powerSeries.getEmbedding(result);
	}

	public int degree() {
		if (degree < 0) {
			if (dimension + ideal.generators().size() == ring.numberOfVariables()) {
				int product = 1;
				for (Polynomial<T> generator : ideal.generators()) {
					product *= generator.degree();
				}
				degree = product;
			} else {
				int accuracy = 0;
				for (Polynomial<T> generator : ideal.generators()) {
					if (generator.degree() > accuracy) {
						accuracy = generator.degree();
					}
				}
				accuracy += 3;
				Rationals q = Rationals.q();
				FormalPowerSeries<Fraction> powerSeries = new FormalPowerSeries<>(q, accuracy);
				PowerSeries<Fraction> oneMinusX = powerSeries.subtract(powerSeries.one(), powerSeries.uniformizer());
				PowerSeries<Fraction> denominator = powerSeries.power(oneMinusX, dimension + 1);
				PowerSeries<Fraction> numerator = powerSeries.multiply(hilbertSeries(powerSeries), denominator);
				PolynomialRing<Fraction> polynomialRing = q.getUnivariatePolynomialRing();
				Fraction result = polynomialRing.evaluate(powerSeries.roundToPolynomial(numerator, accuracy), q.one());
				degree = result.asInteger().intValueExact();
			}
		}
		return degree;
	}

	public boolean isFree() {
		throw new ArithmeticException();
	}

	public List<CoordinateRingElement<T>> getGenerators() {
		List<CoordinateRingElement<T>> result = new ArrayList<>();
		for (int i = 1; i <= this.ring.numberOfVariables(); i++) {
			result.add(getVar(i));
		}
		return result;
	}

	public CoordinateRingElement<T> getVar(int var) {
		return this.getEmbedding(this.ring.getVar(var));
	}

	@Override
	public CoordinateRingElement<T> getEmbedding(Polynomial<T> t) {
		return new CoordinateRingElement<>(this.ideal, ring.getEmbedding(t));
	}

	private CoordinateRingElement<T> getEmbedding(Polynomial<T> t, boolean reduce) {
		return new CoordinateRingElement<>(this.ideal, ring.getEmbedding(t), reduce);
	}

	public CoordinateRingElement<T> getRingEmbedding(T t) {
		return this.getEmbedding(this.ring.getEmbedding(t), false);
	}

	@Override
	public CoordinateRingElement<T> zero() {
		return this.getEmbedding(this.ring.zero(), false);
	}

	@Override
	public CoordinateRingElement<T> one() {
		return this.getEmbedding(this.ring.one(), false);
	}

	@Override
	public CoordinateRingElement<T> add(CoordinateRingElement<T> t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.add(t1.polynomial, t2.polynomial));
	}

	@Override
	public CoordinateRingElement<T> negative(CoordinateRingElement<T> t) {
		return this.getEmbedding(this.ring.negative(t.polynomial), false);
	}

	@Override
	public CoordinateRingElement<T> multiply(CoordinateRingElement<T> t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.multiply(t1.polynomial, t2.polynomial));
	}

	public CoordinateRingElement<T> multiply(T t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.multiply(t1, t2.polynomial));
	}

	public CoordinateRingElement<T> multiplyAssign(T t1, CoordinateRingElement<T> t2) {
		return getEmbedding(this.ring.multiply(t1, t2.polynomial));
	}

	public CoordinateRingElement<T> scalarMultiply(Polynomial<T> t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.multiply(t1, t2.polynomial));
	}

	public CoordinateRingElement<T> scalarMultiplyAssign(Polynomial<T> t1, CoordinateRingElement<T> t2) {
		return getEmbedding(this.ring.multiply(t1, t2.polynomial));
	}

	@Override
	public CoordinateRingElement<T> inverse(CoordinateRingElement<T> t) {
		IdealResult<Polynomial<T>, PolynomialIdeal<T>> asIdeal = ring
				.getIdealWithTransforms(Collections.singletonList(t.getElement()));
		BezoutIdentityResult<Polynomial<T>> bezout = ring.bezoutIdentity(asIdeal.getIdeal(), ideal);
		CoordinateRingElement<T> inverse = zero();
		for (int i = 0; i < asIdeal.getGeneratorExpressions().size(); i++) {
			inverse = add(inverse, getEmbedding(
					ring.multiply(asIdeal.getGeneratorExpressions().get(i).get(0), bezout.getCoeff1().get(i))));
		}
		return inverse;
	}

	@Override
	public boolean isUnit(CoordinateRingElement<T> t) {
		return this.ring.coprime(ring.getIdeal(Collections.singletonList(t.getElement())), ideal);
	}

	public CoordinateRingElement<T> projectToUnit(CoordinateRingElement<T> t) {
		return getEmbedding(ring.projectToUnit(t.polynomial));
	}

	@Override
	public CoordinateRingElement<T> getRandomElement() {
		throw new UnsupportedOperationException();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return BigInteger.valueOf(-1);
	}

	@Override
	public Iterator<CoordinateRingElement<T>> iterator() {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isFinite() {
		return this.ring.getRing().isFinite() && this.krullDimension() == 0;
	}

	@Override
	public BigInteger characteristic() {
		return this.ring.characteristic();
	}

	@Override
	public boolean isIntegral() {
		for (Polynomial<T> generator : ideal.generators()) {
			if (!ring.isIrreducible(generator)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean isReduced() {
		return ideal.isRadical();
	}

	@Override
	public boolean isIrreducible() {
		return ideal.minimalPrimeIdealsOver().size() == 1;
	}

	@Override
	public boolean isZeroDivisor(CoordinateRingElement<T> t) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isCommutative() {
		return ring.isCommutative();
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return false;
	}

	@Override
	public FactorizationResult<CoordinateRingElement<T>, CoordinateRingElement<T>> uniqueFactorization(
			CoordinateRingElement<T> t) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return false;
	}

	@Override
	public boolean isDedekindDomain() {
		return false;
	}

	@Override
	public boolean isDivisible(CoordinateRingElement<T> dividend, CoordinateRingElement<T> divisor) {
		return ring.add(ideal, ring.getIdeal(Collections.singletonList(divisor.getElement())))
				.contains(dividend.getElement());
	}

	@Override
	public QuotientAndRemainderResult<CoordinateRingElement<T>> quotientAndRemainder(CoordinateRingElement<T> dividend,
			CoordinateRingElement<T> divisor) {
		List<Polynomial<T>> generators = new ArrayList<>();
		generators.add(divisor.getElement());
		generators.addAll(ideal.generators());
		IdealResult<Polynomial<T>, PolynomialIdeal<T>> asIdeal = ring.getIdealWithTransforms(generators);
		List<Polynomial<T>> generate = asIdeal.getIdeal().generate(dividend.getElement());
		CoordinateRingElement<T> residue = getEmbedding(asIdeal.getIdeal().residue(dividend.getElement()));
		CoordinateRingElement<T> quotient = zero();
		for (int i = 0; i < generate.size(); i++) {
			quotient = add(quotient,
					getEmbedding(ring.multiply(generate.get(i), asIdeal.getGeneratorExpressions().get(i).get(0))));
		}
		return new QuotientAndRemainderResult<>(quotient, residue);
	}

	@Override
	public BigInteger euclidMeasure(CoordinateRingElement<T> t) {
		return null;
	}

	@Override
	public Iterable<CoordinateRingElement<T>> getUnits() {
		throw new UnsupportedOperationException();
	}

	@Override
	public IdealResult<CoordinateRingElement<T>, CoordinateIdeal<T>> getIdealWithTransforms(
			List<CoordinateRingElement<T>> generators) {
		throw new UnsupportedOperationException();
	}

	public CoordinateIdeal<T> getIdeal(List<CoordinateRingElement<T>> generators) {
		return new CoordinateIdeal<>(this, generators);
	}

	public CoordinateIdeal<T> getIdeal(Ideal<Polynomial<T>> ideal) {
		List<CoordinateRingElement<T>> generators = new ArrayList<>();
		for (Polynomial<T> t : ideal.generators()) {
			generators.add(getEmbedding(t));
		}
		return getIdeal(generators);
	}

	@Override
	public CoordinateIdeal<T> getUnitIdeal() {
		return getIdeal(Collections.singletonList(one()));
	}

	@Override
	public CoordinateIdeal<T> getZeroIdeal() {
		return getIdeal(Collections.emptyList());
	}

	@Override
	public CoordinateIdeal<T> getNilRadical() {
		return getIdeal(ring.radical(getIdeal()));
	}

	@Override
	public CoordinateIdeal<T> intersect(Ideal<CoordinateRingElement<T>> t1, Ideal<CoordinateRingElement<T>> t2) {
		CoordinateIdeal<T> ideal1 = (CoordinateIdeal<T>) t1;
		CoordinateIdeal<T> ideal2 = (CoordinateIdeal<T>) t2;
		PolynomialIdeal<T> intersection = ring.intersect(ideal1.asPolynomialIdeal(), ideal2.asPolynomialIdeal());
		return getIdeal(intersection);
	}

	@Override
	public CoordinateIdeal<T> radical(Ideal<CoordinateRingElement<T>> t) {
		CoordinateIdeal<T> ideal = (CoordinateIdeal<T>) t;
		PolynomialIdeal<T> radical = ring.radical(ideal.asPolynomialIdeal());
		return getIdeal(radical);
	}

	@Override
	public ModuloMaximalIdealResult<CoordinateRingElement<T>, ?> moduloMaximalIdeal(
			Ideal<CoordinateRingElement<T>> ideal) {
		return moduloMaximalIdeal(ideal, ring.moduloMaximalIdeal(((CoordinateIdeal<T>) ideal).asPolynomialIdeal()));
	}

	private <S extends Element<S>> ModuloMaximalIdealResult<CoordinateRingElement<T>, S> moduloMaximalIdeal(
			Ideal<CoordinateRingElement<T>> ideal, ModuloMaximalIdealResult<Polynomial<T>, S> modPolynomialIdeal) {
		return new ModuloMaximalIdealResult<>(this, ideal, modPolynomialIdeal.getField(), new MathMap<>() {
			@Override
			public S evaluate(CoordinateRingElement<T> t) {
				return modPolynomialIdeal.getReduction().evaluate(t.getElement());
			}
		}, new MathMap<>() {
			@Override
			public CoordinateRingElement<T> evaluate(S t) {
				return getEmbedding(modPolynomialIdeal.getLift().evaluate(t));
			}
		});
	}

	@Override
	public ModuloIdealResult<CoordinateRingElement<T>, ?> moduloIdeal(Ideal<CoordinateRingElement<T>> ideal) {
		CoordinateIdeal<T> t = (CoordinateIdeal<T>) ideal;
		return moduloIdeal(t, ring.moduloIdeal(t.asPolynomialIdeal()));
	}

	private <S extends Element<S>> ModuloIdealResult<CoordinateRingElement<T>, S> moduloIdeal(CoordinateIdeal<T> ideal,
			ModuloIdealResult<Polynomial<T>, S> modPolynomialRing) {
		return new ModuloIdealResult<>(this, ideal, modPolynomialRing.getQuotientRing(), new MathMap<>() {
			@Override
			public S evaluate(CoordinateRingElement<T> t) {
				return modPolynomialRing.reduce(t.getElement());
			}
		}, new MathMap<>() {
			@Override
			public CoordinateRingElement<T> evaluate(S t) {
				return getEmbedding(modPolynomialRing.lift(t));
			}
		});
	}

	public CoordinateRingElement<T> substitute(CoordinateRingElement<T> t, List<CoordinateRingElement<T>> values) {
		CoordinateRingElement<T> result = this.zero();
		for (Monomial m : t.polynomial.monomials()) {
			CoordinateRingElement<T> value = this.getRingEmbedding(t.polynomial.coefficient(m));
			int exponents[] = m.exponents();
			for (int i = 0; i < exponents.length; i++) {
				value = this.multiply(value, this.power(values.get(i), exponents[i]));
			}
			result = this.add(result, value);
		}
		return result;
	}

	@Override
	public boolean isGeneratingAlgebra(List<CoordinateRingElement<T>> s) {
		throw new UnsupportedOperationException();
	}

	@Override
	public List<CoordinateRingElement<T>> getAlgebraGenerators() {
		return Collections.singletonList(one());
	}

	@Override
	public boolean isLinearIndependent(List<CoordinateRingElement<T>> s) {
		throw new UnsupportedOperationException();
	}

	@Override
	public List<List<Polynomial<T>>> nonTrivialCombinations(List<CoordinateRingElement<T>> s) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isGeneratingModule(List<CoordinateRingElement<T>> s) {
		throw new UnsupportedOperationException();
	}

	@Override
	public List<CoordinateRingElement<T>> getModuleGenerators() {
		return Collections.singletonList(one());
	}

	@Override
	public Vector<Polynomial<T>> asVector(CoordinateRingElement<T> s) {
		return new Vector<>(Collections.singletonList(s.polynomial));
	}

	public static class CoordinateIdeal<T extends Element<T>> extends AbstractIdeal<CoordinateRingElement<T>> {
		private CoordinateRing<T> ring;
		private List<CoordinateRingElement<T>> generators;
		private SortedMap<Integer, List<Polynomial<T>>> polynomialIdealGeneratorWeights;
		private PolynomialIdeal<T> asPolynomialIdeal;

		private CoordinateIdeal(CoordinateRing<T> ring, List<CoordinateRingElement<T>> generators) {
			super(ring);
			List<Polynomial<T>> generatorPolynomials = new ArrayList<>();
			for (CoordinateRingElement<T> generator : generators) {
				generatorPolynomials.add(generator.polynomial);
			}
			generatorPolynomials.addAll(ring.getIdeal().generators());
			this.ring = ring;
			PolynomialRing<T> p = ring.getPolynomialRing();
			this.asPolynomialIdeal = p.getIdeal(generatorPolynomials);
			this.generators = new ArrayList<>();
			this.polynomialIdealGeneratorWeights = new TreeMap<>();
			for (Polynomial<T> idealGenerator : ring.ideal.generators()) {
				List<Polynomial<T>> generate = asPolynomialIdeal.generate(idealGenerator);
				int index = 0;
				int last = polynomialIdealGeneratorWeights.isEmpty() ? 0
						: polynomialIdealGeneratorWeights.lastKey() + 1;
				boolean found = false;
				for (int i = last; i < generate.size(); i++) {
					if (ring.getIdeal().contains(asPolynomialIdeal.generators().get(i))) {
						continue;
					}
					index = i;
					if (p.isUnit(generate.get(index))) {
						found = true;
						break;
					}
				}
				if (found) {
					Polynomial<T> unit = p.negative(p.inverse(generate.get(index)));
					List<Polynomial<T>> coefficients = new ArrayList<>();
					for (int i = 0; i < generate.size(); i++) {
						if (index == i) {
							coefficients.add(p.zero());
							continue;
						}
						coefficients.add(p.multiply(unit, generate.get(i)));
					}
					polynomialIdealGeneratorWeights.put(index, coefficients);
				}
			}
			for (int i = 0; i < asPolynomialIdeal.generators().size(); i++) {
				if (polynomialIdealGeneratorWeights.containsKey(i)) {
					continue;
				}
				Polynomial<T> generator = asPolynomialIdeal.generators().get(i);
				Polynomial<T> reduced = ring.getIdeal().residue(generator);
				if (!reduced.equals(ring.getPolynomialRing().zero())) {
					this.generators.add(ring.getEmbedding(reduced));
				} else {
					polynomialIdealGeneratorWeights.put(i, Collections.emptyList());
				}
			}
		}

		@Override
		public boolean isPrimary() {
			return asPolynomialIdeal.isPrimary();
		}

		@Override
		public boolean isPrime() {
			return asPolynomialIdeal.isPrime();
		}

		@Override
		public boolean isMaximal() {
			return asPolynomialIdeal.isMaximal();
		}

		public boolean isMaximalOverAlgebraicClosure() {
			return asPolynomialIdeal.isMaximalOverAlgebraicClosure();
		}

		public PolynomialIdeal<T> asPolynomialIdeal() {
			return asPolynomialIdeal;
		}

		@Override
		public List<CoordinateRingElement<T>> generators() {
			return Collections.unmodifiableList(generators);
		}

		@Override
		public List<CoordinateRingElement<T>> generate(CoordinateRingElement<T> t) {
			List<Polynomial<T>> generate = this.asPolynomialIdeal.generate(t.polynomial);
			PolynomialRing<T> p = ring.getPolynomialRing();
			List<Polynomial<T>> preResult = new ArrayList<>();
			for (int i = 0; i < generate.size(); i++) {
				preResult.add(generate.get(i));
			}
			boolean done = false;
			while (!done) {
				done = true;
				for (int i = 0; i < generate.size(); i++) {
					if (polynomialIdealGeneratorWeights.containsKey(i)) {
						List<Polynomial<T>> cofactors = polynomialIdealGeneratorWeights.get(i);
						Polynomial<T> value = preResult.get(i);
						preResult.set(i, p.zero());
						if (cofactors.isEmpty()) {
							continue;
						}
						if (value.equals(p.zero())) {
							continue;
						}
						for (int j = 0; j < generate.size(); j++) {
							if (i == j) {
								continue;
							}
							if (j < i && !cofactors.get(j).equals(p.zero())
									&& polynomialIdealGeneratorWeights.containsKey(j)) {
								done = false;
							}
							preResult.set(j, p.add(preResult.get(j), p.multiply(cofactors.get(j), value)));
						}
					}
				}
			}
			List<CoordinateRingElement<T>> result = new ArrayList<>();
			for (int i = 0; i < preResult.size(); i++) {
				if (!polynomialIdealGeneratorWeights.containsKey(i)) {
					result.add(ring.getEmbedding(preResult.get(i)));
				} else if (!preResult.get(i).equals(p.zero())) {
					throw new ArithmeticException("Did not clear all pre results!");
				}
			}
			return result;
		}

		@Override
		public CoordinateRingElement<T> residue(CoordinateRingElement<T> t) {
			return ring.getEmbedding(asPolynomialIdeal.residue(t.polynomial));
		}

		@Override
		public boolean contains(CoordinateRingElement<T> t) {
			return asPolynomialIdeal.contains(t.polynomial);
		}

		@Override
		public boolean isFinite() {
			return ring.isFinite();
		}

		@Override
		public BigInteger getNumberOfElements() throws InfinityException {
			throw new UnsupportedOperationException();
		}

		public CoordinateRing<T> divideOut() {
			return asPolynomialIdeal().divideOut();
		}

	}
}
