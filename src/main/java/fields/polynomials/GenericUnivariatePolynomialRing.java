package fields.polynomials;

import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.interfaces.Element;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.vectors.FreeModule;
import fields.vectors.Vector;
import util.Pair;
import util.PeekableReader;

public class GenericUnivariatePolynomialRing<T extends Element<T>> extends AbstractPolynomialRing<T> implements
		PolynomialRing<T>/* , DedekindRing<Polynomial<T>, RationalFunction<T>, T> */, UnivariatePolynomialRing<T> {
	private Ring<T> ring;
	private SortedSet<Monomial> monomials;
	private GenericUnivariatePolynomial<T> zero;
	private Map<GenericUnivariatePolynomial<T>, FactorizationResult<Polynomial<T>, T>> squareFreeFactorizatonCache;
	private PolynomialRing<T> asUnivariateBaseRing;

	public GenericUnivariatePolynomialRing(Ring<T> ring) {
		super(false);
		this.ring = ring;
		this.monomials = new TreeSet<>();
		this.monomials.add(getMonomial(new int[] { 0 }));
		this.zero = new GenericUnivariatePolynomial<>(this, Collections.emptyList());
		this.squareFreeFactorizatonCache = new TreeMap<>();
		setVariableName("X");
	}

	@Override
	public UnivariatePolynomialRing<T> withVariableName(String variableName) {
		GenericUnivariatePolynomialRing<T> clone = new GenericUnivariatePolynomialRing<>(ring);
		clone.squareFreeFactorizatonCache = squareFreeFactorizatonCache;
		clone.setVariableName(variableName);
		return clone;
	}

	@Override
	public PolynomialRing<T> withVariableNames(String[] variableNames) {
		return withVariableName(variableNames[0]);
	}

	@Override
	public String getVariableName() {
		return getVariableNames()[0];
	}

	@Override
	public void setVariableName(String variableName) {
		setVariableNames(new String[] { variableName });
	}

	@Override
	public UnivariatePolynomial<T> parse(PeekableReader reader) throws IOException {
		return toUnivariate(super.parse(reader));
	}

	@Override
	public Exactness exactness() {
		return ring.exactness();
	}

	@Override
	public Ring<T> getRing() {
		return ring;
	}

	@Override
	public int numberOfVariables() {
		return 1;
	}

	@Override
	public Comparator<Monomial> getComparator() {
		return Monomial.LEX;
	}

	@Override
	public GenericUnivariatePolynomial<T> zero() {
		return zero;
	}

	@Override
	public GenericUnivariatePolynomial<T> one() {
		return new GenericUnivariatePolynomial<>(this, Collections.singletonList(ring.one()));
	}

	@Override
	public UnivariatePolynomial<T> getVar() {
		return getVarPower(1);
	}

	@Override
	public Polynomial<T> getVar(int i) {
		return getVar();
	}

	@Override
	public UnivariatePolynomial<T> getVarPower(int power) {
		return getEmbedding(ring.one(), power);
	}

	@Override
	public UnivariatePolynomial<T> getVarPower(int i, int power) {
		return getVarPower(power);
	}

	@Override
	public GenericUnivariatePolynomial<T> getEmbedding(T t) {
		return new GenericUnivariatePolynomial<T>(this, Collections.singletonList(t));
	}

	@Override
	public GenericUnivariatePolynomial<T> getPolynomial(List<T> coefficients) {
		return new GenericUnivariatePolynomial<T>(this, coefficients);
	}

	@Override
	@SafeVarargs
	public final UnivariatePolynomial<T> getPolynomial(T... coefficients) {
		return getPolynomial(Arrays.asList(coefficients));
	}

	@Override
	public UnivariatePolynomial<T> getPolynomial(Map<Monomial, T> t) {
		List<T> c = new ArrayList<>();
		for (Monomial m : t.keySet()) {
			while (c.size() <= m.degree()) {
				c.add(ring.zero());
			}
			c.set(m.degree(), t.get(m));
		}
		return getPolynomial(c);
	}

	@Override
	public UnivariatePolynomial<T> getEmbedding(T t, int exponent) {
		List<T> c = new ArrayList<>();
		while (c.size() < exponent) {
			c.add(ring.zero());
		}
		c.add(t);
		return getPolynomial(c);
	}

	@Override
	public UnivariatePolynomial<T> getEmbedding(T t, int[] exponents) {
		return getEmbedding(t, exponents[0]);
	}

	@Override
	public UnivariatePolynomial<T> getLinear(List<T> coeff) {
		if (coeff.size() != 1)
			throw new ArithmeticException("size mismatch!");
		return this.getEmbedding(coeff.get(0), 1);
	}

	@Override
	public GenericUnivariatePolynomial<T> toUnivariate(Polynomial<T> t) {
		if (t instanceof GenericUnivariatePolynomial<?>) {
			return (GenericUnivariatePolynomial<T>) t;
		}
		if (t.numberOfVariables() != 1) {
			throw new ArithmeticException("Number of variables don't match!");
		}
		List<T> c = new ArrayList<>();
		int degree = t.degree();
		for (int i = 0; i <= degree; i++) {
			c.add(t.coefficient(getMonomial(new int[] { i })));
		}
		return getPolynomial(c);
	}

	@Override
	public boolean isHomogeneous(Polynomial<T> t) {
		UnivariatePolynomial<T> p = toUnivariate(t);
		for (int i = 0; i < p.degree(); i++) {
			if (!p.univariateCoefficient(i).equals(ring.zero())) {
				return false;
			}
		}
		return true;
	}

	@Override
	public MultivariatePolynomialRing<T> addVariableWithElimination(int shift) {
		return new MultivariatePolynomialRing<>(ring, shift + 1,
				new Monomial.EliminationOrder(Monomial.LEX, Monomial.LEX, shift));
	}

	@Override
	public UnivariatePolynomial<T> getEmbedding(Polynomial<T> t, int[] map) {
		Map<Monomial, T> c = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			int[] newexp = new int[1];
			for (int i = 0; i < t.numberOfVariables(); i++) {
				if (m.exponents()[i] != 0) {
					newexp[map[i]] = m.exponents()[i];
				}
			}
			c.put(getMonomial(newexp), t.coefficient(m));
		}
		return getPolynomial(c);
	}

	@Override
	public <S extends Element<S>> UnivariatePolynomial<T> getEmbedding(Polynomial<S> t, MathMap<S, T> mapping) {
		UnivariatePolynomial<S> p = ((UnivariatePolynomialRing<S>) t.getPolynomialRing()).toUnivariate(t);
		List<T> c = new ArrayList<>();
		for (int i = 0; i <= p.degree(); i++) {
			c.add(mapping.evaluate(p.univariateCoefficient(i)));
		}
		return getPolynomial(c);
	}

	@Override
	public Polynomial<T> getRandomElement() {
		return getRandomElement(new Random().nextInt(5));
	}

	@Override
	public UnivariatePolynomial<T> getRandomElement(int degree) {
		List<T> c = new ArrayList<>();
		for (int i = 0; i <= degree; i++) {
			c.add(ring.getRandomElement());
		}
		return getPolynomial(c);
	}

	@Override
	public SortedSet<Monomial> monomials(int degree) {
		for (int i = monomials.last().exponents()[0]; i < degree; i++) {
			monomials.add(getMonomial(new int[] { i + 1 }));
		}
		return monomials.headSet(getMonomial(new int[] { degree + 1 }));
	}

	@Override
	public Monomial getMonomial(int[] exponents) {
		return new Monomial(Monomial.LEX, exponents);
	}

	@Override
	public GenericUnivariatePolynomial<T> add(Polynomial<T> s1, Polynomial<T> s2) {
		UnivariatePolynomial<T> p1 = toUnivariate(s1);
		UnivariatePolynomial<T> p2 = toUnivariate(s2);
		int degree = Math.max(s1.degree(), s2.degree());
		List<T> c = new ArrayList<>();
		for (int i = 0; i <= degree; i++) {
			c.add(ring.add(p1.univariateCoefficient(i), p2.univariateCoefficient(i)));
		}
		GenericUnivariatePolynomial<T> result = getPolynomial(c);
		return result;
	}

	@Override
	public UnivariatePolynomial<T> negative(Polynomial<T> s) {
		UnivariatePolynomial<T> p = (UnivariatePolynomial<T>) s;
		List<T> c = new ArrayList<T>();
		int degree = s.degree();
		for (int i = 0; i <= degree; i++) {
			c.add(ring.negative(p.univariateCoefficient(i)));
		}
		return getPolynomial(c);
	}

	@Override
	public GenericUnivariatePolynomial<T> multiply(Polynomial<T> t1, Polynomial<T> t2) {
		GenericUnivariatePolynomial<T> p1 = toUnivariate(t1);
		GenericUnivariatePolynomial<T> p2 = toUnivariate(t2);
		if (p1.degree() < p2.degree()) {
			GenericUnivariatePolynomial<T> tmp = p1;
			p1 = p2;
			p2 = tmp;
		}
		if (p1.degree() <= 0) {
			GenericUnivariatePolynomial<T> result = getPolynomial(
					Collections.singletonList(ring.multiply(p1.univariateCoefficient(0), p2.univariateCoefficient(0))));
			return result;
		}
		int split = p1.degree() / 2 + 1;
		GenericUnivariatePolynomial<T> p1Lsb = getPolynomial(p1.coefficients().subList(0, split));
		GenericUnivariatePolynomial<T> p1Msb = getPolynomial(p1.coefficients().subList(split, p1.degree() + 1));
		if (split > p2.degree()) {
			GenericUnivariatePolynomial<T> lsb = multiply(p1Lsb, p2);
			GenericUnivariatePolynomial<T> msb = multiply(p1Msb, p2);
			GenericUnivariatePolynomial<T> result = add(multiplyPower(split, msb), lsb);
			return result;
		}
		GenericUnivariatePolynomial<T> p2Lsb = getPolynomial(p2.coefficients().subList(0, split));
		GenericUnivariatePolynomial<T> p2Msb = getPolynomial(p2.coefficients().subList(split, p2.degree() + 1));
		UnivariatePolynomial<T> lsb = multiply(p1Lsb, p2Lsb);
		UnivariatePolynomial<T> msb = multiply(p1Msb, p2Msb);
		UnivariatePolynomial<T> p1mixed = add(p1Lsb, p1Msb);
		UnivariatePolynomial<T> p2mixed = add(p2Lsb, p2Msb);
		UnivariatePolynomial<T> mixed = toUnivariate(subtract(multiply(p1mixed, p2mixed), add(lsb, msb)));
		GenericUnivariatePolynomial<T> result = toUnivariate(
				add(lsb, multiplyPower(split, mixed), multiplyPower(2 * split, msb)));
		return result;
	}

	@Override
	public GenericUnivariatePolynomial<T> multiplyPower(int power, Polynomial<T> t) {
		GenericUnivariatePolynomial<T> p = toUnivariate(t);
		List<T> c = new ArrayList<>(power + t.degree() + 1);
		for (int i = 0; i < power; i++) {
			c.add(ring.zero());
		}
		c.addAll(p.coefficients());
		return getPolynomial(c);
	}

	@Override
	public GenericUnivariatePolynomial<T> multiply(T t1, Polynomial<T> t2) {
		return multiplyPower(0, t1, t2);
	}

	@Override
	public GenericUnivariatePolynomial<T> multiplyPower(int power, T t1, Polynomial<T> t2) {
		GenericUnivariatePolynomial<T> p = (GenericUnivariatePolynomial<T>) t2;
		List<T> c = new ArrayList<T>();
		for (int i = 0; i < power; i++) {
			c.add(ring.zero());
		}
		int degree = p.degree();
		for (int i = 0; i <= degree; i++) {
			c.add(ring.multiply(t1, p.univariateCoefficient(i)));
		}
		return getPolynomial(c);
	}

	@Override
	public GenericUnivariatePolynomial<T> divideScalar(Polynomial<T> dividend, T divisor) {
		GenericUnivariatePolynomial<T> p = (GenericUnivariatePolynomial<T>) dividend;
		if (divisor.equals(ring.one())) {
			return p;
		}
		List<T> c = new ArrayList<T>();
		int degree = p.degree();
		for (int i = 0; i <= degree; i++) {
			c.add(ring.divideChecked(p.univariateCoefficient(i), divisor));
		}
		return getPolynomial(c);
	}

	@Override
	public boolean isZeroDivisor(Polynomial<T> t) {
		GenericUnivariatePolynomial<T> p = (GenericUnivariatePolynomial<T>) t;
		if (!ring.isZeroDivisor(p.univariateCoefficient(p.degree()))
				|| !ring.isZeroDivisor(p.univariateCoefficient(0))) {
			return false;
		}
		throw new UnsupportedOperationException("isZeroDivisor not implemented");
	}

	@Override
	public boolean isUnit(Polynomial<T> t) {
		return t.degree() == 0 && ring.isUnit(((UnivariatePolynomial<T>) t).univariateCoefficient(0));
	}

	@Override
	public GenericUnivariatePolynomial<T> inverse(Polynomial<T> t) {
		return getEmbedding(ring.inverse(((UnivariatePolynomial<T>) t).univariateCoefficient(0)));
	}

	@Override
	public boolean isDivisible(Polynomial<T> dividend, Polynomial<T> divisor) {
		if (divisor.equals(zero)) {
			return false;
		}
		return remainder(dividend, divisor).equals(zero());
	}

//	@Override
//	public GeneralQuotientAndRemainderResult<T> generalQuotientAndRemainder(Polynomial<T> dividend,
//			List<Polynomial<T>> divisors) {
//		List<Polynomial<T>> result = new ArrayList<>();
//		GenericUnivariatePolynomial<T> remainder = toUnivariate(dividend);
//		for (Polynomial<T> p : divisors) {
//			GenericUnivariatePolynomial<T> quotient = zero();
//			while (remainder.degree() >= p.degree()
//					&& ring.isDivisible(remainder.leadingCoefficient(), p.leadingCoefficient())) {
//				int exponent = remainder.degree() - p.degree();
//				Ring<T> ringForDivision = ring;
//				if (ring.exactness() == Exactness.FIXED_POINT && ring instanceof LocalField<?, ?>) {
//					@SuppressWarnings("unchecked")
//					LocalField<T, ?> field = (LocalField<T, ?>) ring;
//					int divisorValuation = field.valuation(p.leadingCoefficient()).value();
//					ringForDivision = field.withAccuracy(field.getAccuracy() + Math.abs(divisorValuation));
//				}
//				T factor = ringForDivision.divide(remainder.leadingCoefficient(), p.leadingCoefficient());
//				int deg = remainder.degree();
//				Polynomial<T> approximation = multiplyPower(exponent, factor, p);
//				remainder = toUnivariate(subtract(remainder, approximation));
//				if (remainder.degree() == deg) {
//					throw new ArithmeticException("a * (b / a) != b: " + dividend + "/" + p);
//				}
//				quotient = add(quotient, getEmbedding(factor, exponent));
//			}
//			result.add(quotient);
//		}
//		return new GeneralQuotientAndRemainderResult<>(remainder, result);
//	}

	@Override
	public UnivariatePolynomial<T> pseudoRemainder(Polynomial<T> dividend, Polynomial<T> divisor) {
		T lc = divisor.leadingCoefficient();
		return toUnivariate(
				remainder(multiply(ring.power(lc, dividend.degree() - divisor.degree() + 1), dividend), divisor));
	}

	@Override
	public BigInteger euclidMeasure(Polynomial<T> t) {
		return BigInteger.valueOf(t.degree());
	}

	@Override
	public Iterable<Polynomial<T>> getUnits() {
		return new Iterable<Polynomial<T>>() {

			@Override
			public Iterator<Polynomial<T>> iterator() {
				return new Iterator<Polynomial<T>>() {
					private Iterator<T> it = ring.getUnits().iterator();

					@Override
					public boolean hasNext() {
						return it.hasNext();
					}

					@Override
					public Polynomial<T> next() {
						return getEmbedding(it.next());
					}
				};
			}
		};
	}

	@Override
	public FactorizationResult<Polynomial<T>, T> squareFreeFactorization(Polynomial<T> t) {
		GenericUnivariatePolynomial<T> f = (GenericUnivariatePolynomial<T>) t;
		if (this.squareFreeFactorizatonCache.containsKey(f)) {
			return this.squareFreeFactorizatonCache.get(f);
		}
		SortedMap<Polynomial<T>, Integer> result = new TreeMap<>();
		GenericUnivariatePolynomial<T> tPrime = derivative(t);
		// Over algebraic closure, this decomposes into the factors of t with one degree
		// less.
		GenericUnivariatePolynomial<T> degreeMinusI = gcd(t, tPrime);
		// Over the algebraic closure, this decomposes into the linear factors of t.
		GenericUnivariatePolynomial<T> linearFactors = toUnivariate(this.divideChecked(t, degreeMinusI));
		int i = 0;
		// Finding all factors of degree i with gcd(i, char F) = 1.
		T unit = ring.one();
		while (linearFactors.degree() != 0) {
			i++;
			// Over the algebraic closure, this is all the factors with multiplicity higher
			// than i.
			GenericUnivariatePolynomial<T> linearFactorsHigherDegree = this.gcd(linearFactors, degreeMinusI);
			// This are all the linear factors that have multiplicity i in t.
			GenericUnivariatePolynomial<T> factorThisDegree = toUnivariate(
					this.divideChecked(linearFactors, linearFactorsHigherDegree));
			if (factorThisDegree.degree() > 0) {
				result.put(factorThisDegree, i);
			} else {
				unit = ring.multiply(factorThisDegree.leadingCoefficient(), unit);
			}
			// Updating the linear factors.
			linearFactors = linearFactorsHigherDegree;
			// Removing multiplicity i factors.
			degreeMinusI = toUnivariate(this.divideChecked(degreeMinusI, linearFactorsHigherDegree));
		}
		// Remaining factors are multiplicity kp
		if (degreeMinusI.degree() != 0) {
			if (hasCharacteristicRoot(degreeMinusI)) {
				Polynomial<T> pthRoot = characteristicRoot(degreeMinusI);
				FactorizationResult<Polynomial<T>, T> factors = this.squareFreeFactorization(pthRoot);
				for (Polynomial<T> factor : factors.primeFactors()) {
					int power = this.characteristic().intValueExact() * factors.multiplicity(factor);
					degreeMinusI = toUnivariate(divideChecked(degreeMinusI, power(factor, power)));
					result.put(factor, power);
				}
			} else {
				throw new ArithmeticException("Characteristic zero should not have any remaining factors");
			}
		}
		if (!isUnit(degreeMinusI)) {
			throw new ArithmeticException("Polynomial was not content free!");
		}
		FactorizationResult<Polynomial<T>, T> factors = new FactorizationResult<>(
				ring.multiply(unit, degreeMinusI.leadingCoefficient(), linearFactors.leadingCoefficient()), result);
		this.squareFreeFactorizatonCache.put(f, factors);
		return factors;
	}

	@Override
	public List<Polynomial<T>> getAlgebraGenerators() {
		return Collections.singletonList(getVar());
	}

	@Override
	public Vector<T> asVector(Polynomial<T> s) {
		return asVector(s, s.degree());
	}

	@Override
	public Vector<T> asVector(Polynomial<T> s, int degree) {
		if (s.degree() > degree) {
			throw new ArithmeticException("degree too low");
		}
		UnivariatePolynomial<T> p = (UnivariatePolynomial<T>) s;
		List<T> vector = new ArrayList<T>();
		for (int i = 0; i <= degree; i++) {
			vector.add(p.univariateCoefficient(i));
		}
		return new Vector<T>(vector);
	}

	@Override
	public GenericUnivariatePolynomial<T> fromVector(Vector<T> vector) {
		return getPolynomial(vector.asList());
	}

	@Override
	public Vector<T> asVector(Polynomial<T> t, int[] degrees) {
		return asVector(t, degrees[0]);
	}

	@Override
	public Polynomial<T> fromVector(Vector<T> t, int[] degrees) {
		return fromVector(t);
	}

	@Override
	public boolean isLinearIndependent(List<Polynomial<T>> s) {
		int degree = 0;
		for (Polynomial<T> p : s) {
			degree = Math.max(degree, p.degree());
		}
		FreeModule<T> free = new FreeModule<>(getRing(), degree + 1);
		List<Vector<T>> asVectors = new ArrayList<>();
		for (Polynomial<T> p : s) {
			asVectors.add(asVector(p, degree));
		}
		return free.isLinearIndependent(asVectors);
	}

	@Override
	public List<Vector<T>> nonTrivialCombinations(List<Polynomial<T>> s) {
		int degree = 0;
		for (Polynomial<T> p : s) {
			degree = Math.max(degree, p.degree());
		}
		FreeModule<T> free = new FreeModule<>(getRing(), degree + 1);
		List<Vector<T>> asVectors = new ArrayList<>();
		for (Polynomial<T> p : s) {
			asVectors.add(asVector(p, degree));
		}
		return free.nonTrivialCombinations(asVectors);
	}

	@Override
	public Iterator<Polynomial<T>> iterator() {
		return new Iterator<Polynomial<T>>() {
			private int degree = -1;
			private Iterator<UnivariatePolynomial<T>> it = null;

			@Override
			public boolean hasNext() {
				return true;
			}

			@Override
			public Polynomial<T> next() {
				if (degree < 0) {
					degree = 0;
					it = polynomials(degree);
					return zero();
				}
				if (!it.hasNext()) {
					degree++;
					it = polynomials(degree);
				}
				return it.next();
			}
		};
	}

	@Override
	public Iterator<UnivariatePolynomial<T>> polynomials(int degree) {
		if (degree < 0) {
			return new Iterator<UnivariatePolynomial<T>>() {
				private boolean used = false;

				@Override
				public boolean hasNext() {
					return !used;
				}

				@Override
				public GenericUnivariatePolynomial<T> next() {
					if (used) {
						throw new RuntimeException("Iterator at end!");
					}
					used = true;
					return zero();
				}
			};
		}
		return new Iterator<UnivariatePolynomial<T>>() {
			private int limitIndex;
			private long limit;
			private long[] index;
			private List<T> current = new ArrayList<>();
			private List<Iterator<T>> its = null;

			@Override
			public boolean hasNext() {
				return current != null;
			}

			@Override
			public GenericUnivariatePolynomial<T> next() {
				if (its == null) {
					its = new ArrayList<>();
					limit = 0;
					index = new long[degree + 1];
					limitIndex = degree;
					for (int i = 0; i <= degree; i++) {
						its.add(ring.iterator());
						current.add(its.get(i).next());
						index[i] = 0;
					}
				}
				if (current == null) {
					throw new RuntimeException("Iterator at end!");
				}
				List<T> copy = new ArrayList<>();
				copy.addAll(current);
				GenericUnivariatePolynomial<T> next = getPolynomial(copy);
				boolean broken = false;
				boolean hasNext = false;
				for (int i = 0; i <= degree; i++) {
					if (its.get(i).hasNext()) {
						hasNext = true;
						break;
					}
				}
				while (!broken && hasNext) {
					for (int i = 0; i <= degree; i++) {
						if (!its.get(i).hasNext()) {
							continue;
						}
						if (i == limitIndex) {
							continue;
						}
						if (i < limitIndex && index[i] >= limit) {
							continue;
						}
						if (i > limitIndex && index[i] >= limit - 1) {
							continue;
						}
						current.set(i, its.get(i).next());
						index[i]++;
						for (int j = 0; j < i; j++) {
							if (j == limitIndex) {
								continue;
							}
							its.set(j, ring.iterator());
							current.set(j, its.get(j).next());
							index[j] = 0;
						}
						broken = true;
						break;
					}
					if (!broken) {
						if (limitIndex < degree) {
							limitIndex++;
							current.set(limitIndex, its.get(limitIndex).next());
							index[limitIndex]++;
							for (int i = 0; i <= degree; i++) {
								if (i == limitIndex) {
									continue;
								}
								its.set(i, ring.iterator());
								current.set(i, its.get(i).next());
								index[i] = 0;

							}
						} else {
							limit++;
							limitIndex = 0;
							current.set(0, its.get(0).next());
							index[0]++;
							for (int i = 1; i <= degree; i++) {
								its.set(i, ring.iterator());
								current.set(i, its.get(i).next());
								index[i] = 0;
							}
						}
						break;
					}
				}
				if (!hasNext) {
					current = null;
				}
				return next;
			}
		};
	}

	@Override
	public Iterable<UnivariatePolynomial<T>> polynomialSet(int degree) {
		return new Iterable<UnivariatePolynomial<T>>() {

			@Override
			public Iterator<UnivariatePolynomial<T>> iterator() {
				return polynomials(degree);
			}
		};
	}

	@Override
	public Iterator<UnivariatePolynomial<T>> monicPolynomials(int degree) {
		return new Iterator<UnivariatePolynomial<T>>() {
			private Iterator<UnivariatePolynomial<T>> it = polynomials(degree - 1);

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public GenericUnivariatePolynomial<T> next() {
				return add(it.next(), getVarPower(degree));
			}
		};
	}

	@Override
	public Iterable<UnivariatePolynomial<T>> monicPolynomialSet(int degree) {
		return new Iterable<UnivariatePolynomial<T>>() {

			@Override
			public Iterator<UnivariatePolynomial<T>> iterator() {
				return monicPolynomials(degree);
			}
		};
	}
	
	@Override
	public ChineseRemainderPreparation<Polynomial<T>> prepareInterpolation(List<T> interpolationPoints) {
		List<Polynomial<T>> linear = new ArrayList<>();
		for (T point : interpolationPoints) {
			linear.add(subtract(getVar(), getEmbedding(point)));
		}
		return prepareChineseRemainderTheoremModuli(linear);
	}
	
	@Override
	public UnivariatePolynomial<T> interpolate(ChineseRemainderPreparation<Polynomial<T>> preparation, List<T> interpolationValues) {
		List<Polynomial<T>> constant = new ArrayList<>();
		for (T point : interpolationValues) {
			constant.add(getEmbedding(point));
		}
		return toUnivariate(chineseRemainderTheorem(constant, preparation));
	}
	
	@Override
	public UnivariatePolynomial<T> interpolate(List<T> interpolationPoints, List<T> interpolationValues) {
		return interpolate(prepareInterpolation(interpolationPoints), interpolationValues);
	}

	@Override
	public final T evaluate(Polynomial<T> t, List<T> ts) {
		if (ts.size() != 1) {
			throw new ArithmeticException("Only a single variable needed");
		}
		GenericUnivariatePolynomial<T> p = (GenericUnivariatePolynomial<T>) t;
		if (p.degree() < 0) {
			return ring.zero();
		}
		T value = ts.get(0);
		T result = p.univariateCoefficient(p.degree());
		for (int i = p.degree() - 1; i >= 0; i--) {
			result = ring.multiply(result, value);
			result = ring.add(result, p.univariateCoefficient(i));
		}
		return result;
	}

	@Override
	public Polynomial<T> partiallyEvaluate(Polynomial<T> t, List<T> ts) {
		if (ts.get(0) == null) {
			return t;
		}
		return getEmbedding(evaluate(t, ts));
	}

	@Override
	public GenericUnivariatePolynomial<T> substitute(Polynomial<T> t, List<Polynomial<T>> values) {
		GenericUnivariatePolynomial<T> p = toUnivariate(t);
		if (values.get(0).equals(getVar())) {
			return p;
		}
		GenericUnivariatePolynomial<T> value = toUnivariate(values.get(0));
		GenericUnivariatePolynomial<T> result = getEmbedding(p.univariateCoefficient(p.degree()));
		for (int i = p.degree() - 1; i >= 0; i--) {
			result = multiply(result, value);
			result = add(result, getEmbedding(p.univariateCoefficient(i)));
		}
		return result;
	}

	@Override
	public boolean canBeNormalized(Polynomial<T> t) {
		return ring.isDivisible(content(t), t.leadingCoefficient());
	}

	@Override
	public GenericUnivariatePolynomial<T> normalize(Polynomial<T> t) {
		return divideScalar(t, t.leadingCoefficient());
	}

	@Override
	public T content(Polynomial<T> t) {
		GenericUnivariatePolynomial<T> p = toUnivariate(t);
		if (p.content() != null) {
			return p.content();
		}
		T content = t.leadingCoefficient();
		for (int i = 0; i < p.degree(); i++) {
			content = ring.gcd(content, p.univariateCoefficient(i));
		}
		p.setContent(content);
		return content;
	}

	@Override
	public GenericUnivariatePolynomial<T> contentFree(Polynomial<T> t) {
		return divideScalar(t, content(t));
	}

	@Override
	public GenericUnivariatePolynomial<T> depress(Polynomial<T> t) {
		GenericUnivariatePolynomial<T> polynomial = toUnivariate(t);
		T coefficient = polynomial.univariateCoefficient(polynomial.degree() - 1);
		T degree = ring.getInteger(polynomial.degree());
		if (coefficient.equals(ring.zero())) {
			return polynomial;
		}
		if (!ring.isDivisible(coefficient, degree)) {
			return polynomial;
		}
		Polynomial<T> substitute = subtract(getVar(), getEmbedding(ring.divideChecked(coefficient, degree)));
		return substitute(polynomial, Collections.singletonList(substitute));
	}

	@Override
	public GenericUnivariatePolynomial<T> derivative(Polynomial<T> t, int variable) {
		return derivative(t);
	}

	@Override
	public GenericUnivariatePolynomial<T> derivative(Polynomial<T> t) {
		GenericUnivariatePolynomial<T> p = (GenericUnivariatePolynomial<T>) t;
		List<T> c = new ArrayList<>();
		for (int i = 0; i < p.degree(); i++) {
			c.add(ring.multiply(i + 1, p.univariateCoefficient(i + 1)));
		}
		return getPolynomial(c);
	}

	@Override
	public Vector<T> trace(Polynomial<T> t, int degree) {
		return toUnivariate(t).trace(degree);
	}

	@Override
	public T resultant(Polynomial<T> t1, Polynomial<T> t2) {
		return extendedResultant(t1, t2).getResultant();
	}

	@Override
	public Polynomial<T> resultant(Polynomial<T> t1, Polynomial<T> t2, int variable) {
		if (variable != 1) {
			throw new ArithmeticException("wrong number of variables");
		}
		return eliminateVariable().getEmbedding(resultant(t1, t2));
	}

	@Override
	public PolynomialRing<T> eliminateVariable() {
		if (asUnivariateBaseRing != null) {
			return asUnivariateBaseRing;
		}
		asUnivariateBaseRing = AbstractPolynomialRing.getPolynomialRing(ring, 0, Monomial.GREVLEX);
		return asUnivariateBaseRing;
	}

	@Override
	public UnivariatePolynomial<Polynomial<T>> asUnivariatePolynomial(Polynomial<T> t, int variable) {
		if (variable != 1) {
			throw new ArithmeticException("wrong number of variables");
		}
		return eliminateVariable().getUnivariatePolynomialRing().getEmbedding(t, new MathMap<>() {
			@Override
			public Polynomial<T> evaluate(T t) {
				return eliminateVariable().getEmbedding(t);
			}
		});
	}

	@Override
	public Polynomial<T> fromUnivariatePolynomial(UnivariatePolynomial<Polynomial<T>> t, int variable) {
		return getEmbedding(t, new MathMap<>() {
			@Override
			public T evaluate(Polynomial<T> t) {
				return t.coefficient(eliminateVariable().getMonomial(new int[] {}));
			}
		});
	}

	@SuppressWarnings("unchecked")
	@Override
	public GenericUnivariatePolynomial<T> gcd(Polynomial<T> t1, Polynomial<T> t2) {
		if (t1.equals(zero())) {
			return toUnivariate(t2);
		}
		if (t2.equals(zero())) {
			return toUnivariate(t1);
		}
		if (ring.exactness() == Exactness.FIXED_POINT && ring instanceof DiscreteValuationRing<?, ?>) {
			return gcdLocalRing(t1, t2, (DiscreteValuationRing<T, ?>) ring);
		}
		if (ring.exactness() == Exactness.FIXED_POINT && ring instanceof DiscreteValuationField<?, ?>) {
			return toUnivariate(super.gcd(t1, t2));
		}
		T content1 = content(t1);
		T content2 = content(t2);
		T content = ring.gcd(content1, content2);
		Polynomial<T> contentFree = contentFree(extendedResultant(contentFree(t1), contentFree(t2)).getGcd());
		return multiply(content, contentFree);
	}

	private <S extends Element<S>> GenericUnivariatePolynomial<T> gcdLocalRing(Polynomial<T> t1, Polynomial<T> t2,
			DiscreteValuationRing<T, S> ring) {
		T content1 = content(t1);
		T content2 = content(t2);
		T content = ring.gcd(content1, content2);
		Polynomial<T> contentFree = contentFree(
				ring.localField().getUnivariatePolynomialRing().gcd(contentFree(t1), contentFree(t2)));
		return multiply(content, contentFree);
	}

	@Override
	public ExtendedEuclideanResult<Polynomial<T>> extendedEuclidean(Polynomial<T> t1, Polynomial<T> t2) {
		if (!isEuclidean()) {
			throw new ArithmeticException("Ring " + ring + " is not Euclidean!");
		}
		if (ring.exactness() == Exactness.FIXED_POINT) {
			return super.extendedEuclidean(t1, t2);
		}
		ExtendedResultantResult<T> er = extendedResultant(t1, t2);
		return new ExtendedEuclideanResult<>(er.getGcd(), er.getCoeff1(), er.getCoeff2());
	}

	@Override
	public T discriminant(Polynomial<T> t1) {
		return resultant(t1, derivative(t1));
	}

	@Override
	public ExtendedResultantResult<T> extendedResultant(Polynomial<T> t1, Polynomial<T> t2) {
		if (t1.degree() < t2.degree()) {
			ExtendedResultantResult<T> swapped = extendedResultant(t2, t1);
			ExtendedResultantResult<T> result = new ExtendedResultantResult<>(swapped.getResultant(),
					swapped.getResultantCoeff2(), swapped.getResultantCoeff1(), swapped.getGcd(), swapped.getCoeff2(),
					swapped.getCoeff1());
			return result;
		}
		if (t2.degree() == 0) {
			if (t1.degree() == 0) {
				Pair<T, T> bezout = ring.bezoutIdentity(t1.leadingCoefficient(), t2.leadingCoefficient());
				return new ExtendedResultantResult<>(ring.one(), getEmbedding(bezout.getFirst()),
						getEmbedding(bezout.getSecond()), one(), getEmbedding(bezout.getFirst()),
						getEmbedding(bezout.getSecond()));
			}
			T power = ring.power(t2.leadingCoefficient(), t1.degree() - 1);
			T resultant = ring.multiply(power, t2.leadingCoefficient());
			ExtendedResultantResult<T> result = new ExtendedResultantResult<>(resultant, zero(), getEmbedding(power),
					getEmbedding(resultant), zero(), getEmbedding(power));
			return result;
		}
		Polynomial<T> rPrevPrev = null;
		Polynomial<T> rPrev = t1;
		Polynomial<T> coeff1Prev = one();
		Polynomial<T> coeff2Prev = zero();
		Polynomial<T> r = t2;
		Polynomial<T> coeff1 = zero();
		Polynomial<T> coeff2 = one();
		T gamma = t1.leadingCoefficient();
		int degree = 0;
		T psi = null;
		T beta;
		for (int i = 1; !r.equals(zero()); i++) {
			T gammaPrev = gamma;
			gamma = r.leadingCoefficient();
			int degreePrev = degree;
			degree = rPrev.degree() - r.degree();
			if (i == 1) {
				beta = ring.getInteger(degree % 2 == 0 ? -1 : 1);
				psi = ring.getInteger(-1);
			} else {
				T psiPrev = psi;
				psi = ring.divide(ring.power(ring.negative(gammaPrev), degreePrev),
						ring.power(psiPrev, degreePrev - 1));
				beta = ring.negative(ring.multiply(gammaPrev, ring.power(psi, degree)));
			}
			T gammaPower = ring.power(gamma, degree + 1);
			QuotientAndRemainderResult<Polynomial<T>> qr = quotientAndRemainder(multiply(gammaPower, rPrev), r);
			Polynomial<T> rNext = divideScalar(qr.getRemainder(), beta);
			// rPrev = coeff1Prev * t1 + coeff2Prev * t2
			// r = coeff1 * t1 + coeff2 * t2
			// rNext = rem(gamma^(degree + 1) rPrev, r) / beta
			// gamma^(degree + 1) rPrev = q * r + rNext * beta
			// rNext * beta = gamma^(degree + 1) rPrev - q * r
			// rNext * beta = gamma^(degree + 1) (c1Prev * t1 + c2Prev * t2) - q * (c1 * t1
			// + c2 * t2)
			// rNext * beta = (gamma^(degree + 1) c1Prev - q * c1) * t1 + (gamma^(degree +
			// 1) c2Prev - q * c2) * t2
			// rNext = (gamma^(degree + 1) c1Prev - q * c1) / beta * t1 + (gamma^(degree +
			// 1) c2Prev - q * c2) / beta * t2
			Polynomial<T> coeff1Next = divideScalar(
					subtract(multiply(gammaPower, coeff1Prev), multiply(qr.getQuotient(), coeff1)), beta);
			Polynomial<T> coeff2Next = divideScalar(
					subtract(multiply(gammaPower, coeff2Prev), multiply(qr.getQuotient(), coeff2)), beta);
			rPrevPrev = rPrev;
			rPrev = r;
			coeff1Prev = coeff1;
			coeff2Prev = coeff2;
			r = rNext;
			coeff1 = coeff1Next;
			coeff2 = coeff2Next;
		}
		UnivariatePolynomial<T> gcd = toUnivariate(rPrev);
		T resultant;
		UnivariatePolynomial<T> resultantCoeff1;
		UnivariatePolynomial<T> resultantCoeff2;
		if (rPrevPrev != null && rPrev.degree() == 0) {
			T powerM1 = ring.power(rPrev.leadingCoefficient(), rPrevPrev.degree() - 1);
			resultant = ring.multiply(powerM1, rPrev.leadingCoefficient());
			resultantCoeff1 = multiply(powerM1, coeff1Prev);
			resultantCoeff2 = multiply(powerM1, coeff2Prev);
		} else {
			resultant = ring.zero();
			resultantCoeff1 = toUnivariate(t2);
			resultantCoeff2 = negative(t1);
		}
		ExtendedResultantResult<T> result = new ExtendedResultantResult<>(resultant, resultantCoeff1, resultantCoeff2,
				gcd, toUnivariate(coeff1Prev), toUnivariate(coeff2Prev));
		return result;
	}

	@Override
	public GenericUnivariatePolynomial<T> round(Polynomial<T> t, int degree) {
		GenericUnivariatePolynomial<T> p = toUnivariate(t);
		return getPolynomial(p.coefficients().subList(0, Math.min(degree, p.degree() + 1)));
	}
}
