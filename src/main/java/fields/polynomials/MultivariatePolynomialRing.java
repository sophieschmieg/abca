package fields.polynomials;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.polynomials.Monomial.EliminationOrder;
import fields.vectors.FreeModule;
import fields.vectors.Vector;

public class MultivariatePolynomialRing<T extends Element<T>> extends AbstractPolynomialRing<T>
		implements PolynomialRing<T> {
	private Ring<T> ring;
	private int numvars;
	private Comparator<Monomial> comparator;
	private PolynomialRing<T> eliminatedVariable;

	MultivariatePolynomialRing(Ring<T> ring, int numvars, Comparator<Monomial> comparator) {
		super(true);
		this.ring = ring;
		this.numvars = numvars;
		this.comparator = comparator;
		setVariableNames(getMonomial(new int[numberOfVariables()]).defaultVariableNames());
	}

	@Override
	public Exactness exactness() {
		return ring.exactness();
	}

	public Ring<T> getRing() {
		return this.ring;
	}

	@Override
	public PolynomialRing<T> withVariableNames(String[] variableNames) {
		PolynomialRing<T> clone = new MultivariatePolynomialRing<>(ring, numvars, comparator);
		clone.setVariableNames(variableNames);
		return clone;
	}

	@Override
	public boolean isLinearIndependent(List<Polynomial<T>> s) {
		Monomial degrees = getMonomial(new int[numvars]);
		for (Polynomial<T> t : s) {
			degrees = degrees.lcm(t.leadingMonomial());
		}
		int dimension = 1;
		int[] dimensions = new int[numvars];
		for (int i = 0; i < numvars; i++) {
			dimensions[i] = degrees.exponents()[i] + 1;
			dimension *= dimensions[i];
		}
		FreeModule<T> free = new FreeModule<>(getRing(), dimension);
		List<Vector<T>> asVectors = new ArrayList<>();
		for (Polynomial<T> t : s) {
			asVectors.add(asVector((MultivariatePolynomial<T>) t, dimensions));
		}
		return free.isLinearIndependent(asVectors);
	}

	@Override
	public List<Vector<T>> nonTrivialCombinations(List<Polynomial<T>> s) {
		Monomial degrees = getMonomial(new int[numvars]);
		for (Polynomial<T> t : s) {
			degrees = degrees.lcm(t.leadingMonomial());
		}
		int dimension = 1;
		int[] dimensions = new int[numvars];
		for (int i = 0; i < numvars; i++) {
			dimensions[i] = degrees.exponents()[i] + 1;
			dimension *= dimensions[i];
		}
		FreeModule<T> free = new FreeModule<>(getRing(), dimension);
		List<Vector<T>> asVectors = new ArrayList<>();
		for (Polynomial<T> t : s) {
			asVectors.add(asVector((MultivariatePolynomial<T>) t, dimensions));
		}
		return free.nonTrivialCombinations(asVectors);
	}

	public MultivariatePolynomial<T> toMultivariate(Polynomial<T> t) {
		if (t instanceof MultivariatePolynomial<?> && t.getPolynomialRing() == this) {
			return (MultivariatePolynomial<T>) t;
		}
		Map<Monomial, T> c = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			c.put(m, t.coefficient(m));
		}
		return getPolynomial(c);
	}

	@Override
	public Vector<T> asVector(Polynomial<T> s, int[] dimensions) {
		int dimension = 1;
		for (int i = 0; i < numvars; i++) {
			if (dimensions[i] < s.degree(i + 1) + 1) {
				throw new ArithmeticException("Not high dimensional enough");
			}
			dimension *= dimensions[i];
		}
		List<T> asList = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			int value = i;
			int[] exponents = new int[numvars];
			for (int j = 0; j < numvars; j++) {
				exponents[j] = value % dimensions[j];
				value /= dimensions[j];
			}
			asList.add(s.coefficient(getMonomial(exponents)));
		}
		return new Vector<>(asList);
	}

	@Override
	public Vector<T> asVector(Polynomial<T> s) {
		int[] dimensions = new int[numvars];
		for (int i = 0; i < numvars; i++) {
			dimensions[i] = s.degree(i + 1) + 1;
		}
		return asVector(s, dimensions);
	}

	@Override
	public Polynomial<T> fromVector(Vector<T> t, int[] degrees) {
		Map<Monomial, T> result = new TreeMap<>();
		for (int i = 0; i < t.dimension(); i++) {
			int value = i;
			int[] exponents = new int[numvars];
			for (int j = 0; j < numvars; j++) {
				exponents[j] = value % degrees[j];
				value /= degrees[j];
			}
			result.put(getMonomial(exponents), t.get(i + 1));
		}
		return getPolynomial(result);
	}

	@Override
	public List<Polynomial<T>> getAlgebraGenerators() {
		List<Polynomial<T>> result = new ArrayList<>();
		for (int i = 1; i <= numvars; i++) {
			result.add(getVar(i));
		}
		return result;
	}

	public MultivariatePolynomial<T> getEmbedding(T t, int[] exponents) {
		return this.getEmbedding(t, getMonomial(exponents));
	}

	public MultivariatePolynomial<T> getEmbedding(T t, Monomial monomial) {
		return getPolynomial(Collections.singletonMap(monomial, t));
	}

	public MultivariatePolynomial<T> getEmbedding(T t) {
		return this.getEmbedding(t, new int[numvars]);
	}

	@Override
	public MultivariatePolynomial<T> getEmbedding(Polynomial<T> t, int[] map) {
		Map<Monomial, T> c = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			int[] newexp = new int[this.numvars];
			for (int i = 0; i < t.numberOfVariables(); i++) {
				if (m.exponents()[i] != 0)
					newexp[map[i]] = m.exponents()[i];
			}
			c.put(getMonomial(newexp), t.coefficient(m));
		}
		return getPolynomial(c);
	}

	@Override
	public <S extends Element<S>> Polynomial<T> getEmbedding(Polynomial<S> t, MathMap<S, T> mapping) {
		Map<Monomial, T> coeff = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			coeff.put(m, mapping.evaluate(t.coefficient(m)));
		}
		return getPolynomial(coeff);
	}

	@Override
	public MultivariatePolynomial<T> getPolynomial(Map<Monomial, T> coeff) {
		return new MultivariatePolynomial<T>(coeff, this);
	}

	@SafeVarargs
	public final MultivariatePolynomial<T> getUnivariatePolynomial(T... coeffs) {
		Map<Monomial, T> coeffMap = new TreeMap<>(comparator);
		for (int i = 0; i < coeffs.length; i++) {
			coeffMap.put(getMonomial(new int[] { coeffs.length - 1 - i }), coeffs[i]);
		}
		return getPolynomial(coeffMap);
	}

	public MultivariatePolynomial<T> getVar(int var) {
		return this.getVarPower(var, 1);
	}

	public MultivariatePolynomial<T> getVarPower(int var, int power) {
		int[] exponents = new int[numvars];
		if (var < 1 || var > numvars)
			throw new ArithmeticException("wrong number of vars");
		var--;
		exponents[var] = power;
		return this.getEmbedding(this.ring.one(), exponents);
	}

	@Override
	public MultivariatePolynomial<T> getLinear(List<T> coeff) {
		if (this.numvars != coeff.size())
			throw new ArithmeticException("size mismatch!");
		Map<Monomial, T> coeffs = new TreeMap<>();
		int[] exponents = new int[numvars];
		for (int i = 0; i < numvars; i++) {
			exponents[i] = 1;
			if (i > 0)
				exponents[i - 1] = 0;
			coeffs.put(getMonomial(exponents), coeff.get(i));
		}
		return this.getPolynomial(coeffs);
	}

	public PolynomialRing<T> eliminateVariable() {
		if (eliminatedVariable != null) {
			return eliminatedVariable;
		}
		eliminatedVariable = AbstractPolynomialRing.getPolynomialRing(ring, numvars - 1, Monomial.GREVLEX);
		return eliminatedVariable;
	}

	public UnivariatePolynomial<Polynomial<T>> asUnivariatePolynomial(Polynomial<T> t, int variable) {
		if (variable < 1 || variable > numvars) {
			throw new ArithmeticException("Variable out of bounds");
		}
		PolynomialRing<T> eliminated = eliminateVariable();
		UnivariatePolynomialRing<Polynomial<T>> eliminatedRing = eliminated.getUnivariatePolynomialRing();
		UnivariatePolynomial<Polynomial<T>> result = eliminatedRing.zero();
		for (Monomial m : t.monomials()) {
			int[] exponents = new int[numvars - 1];
			for (int i = 0; i < numvars; i++) {
				if (i == variable - 1) {
					continue;
				}
				exponents[i < variable - 1 ? i : i - 1] = m.exponents()[i];
			}
			result = eliminatedRing.add(result,
					eliminatedRing.multiply(eliminated.getEmbedding(t.coefficient(m), exponents),
							eliminatedRing.getVarPower(m.exponents()[variable - 1])));
		}
		return result;
	}

	public Polynomial<T> fromUnivariatePolynomial(UnivariatePolynomial<Polynomial<T>> t, int variable) {
		if (variable < 1 || variable > numvars) {
			throw new ArithmeticException("Variable out of bounds");
		}
		int[] map = new int[numvars];
		for (int i = 0; i < numvars - 1; i++) {
			map[i] = i < variable - 1 ? i : i + 1;
		}
		Polynomial<T> result = zero();
		for (int i = 0; i <= t.degree(); i++) {
			result = add(result, multiply(getEmbedding(t.univariateCoefficient(i), map), getVarPower(variable, i)));
		}
		return result;
	}

	public MultivariatePolynomialRing<T> addVariableWithElimination(int shift) {
		if (shift < 0)
			throw new ArithmeticException("Cannot count");
		return new MultivariatePolynomialRing<T>(ring, this.numvars + shift,
				new EliminationOrder(Monomial.LEX, this.comparator, shift));
	}

	@Override
	public MultivariatePolynomial<T> getRandomElement() {
		int degree = (int) Math.abs(new Random().nextGaussian() * 2.0);
		return getRandomElement(degree);
	}

	private MultivariatePolynomial<T> getRandomElement(int degree, int var) {
		if (degree < 0) {
			return this.zero();
		}
		if (var > numvars) {
			return this.getEmbedding(ring.getRandomElement());
		}
		MultivariatePolynomial<T> result = this.zero();
		for (int i = 0; i <= degree; i++) {
			result = this.add(result,
					this.multiply(this.getVarPower(var, i), this.getRandomElement(degree - i, var + 1)));
		}
		return result;
	}

	@Override
	public MultivariatePolynomial<T> getRandomElement(int degree) {
		return this.getRandomElement(degree, 1);
	}

	@Override
	public Iterator<Polynomial<T>> iterator() {
		throw new InfinityException();
	}

	@Override
	public MultivariatePolynomial<T> zero() {
		return this.getEmbedding(this.ring.zero());
	}

	@Override
	public MultivariatePolynomial<T> add(Polynomial<T> t1, Polynomial<T> t2) {
		Map<Monomial, T> coeff = new TreeMap<>();
		coeff.put(getMonomial(new int[this.numvars]), this.ring.zero());
		Set<Monomial> monomeset = new TreeSet<>();
		monomeset.addAll(t1.monomials());
		monomeset.addAll(t2.monomials());
		for (Monomial m : monomeset) {
			coeff.put(m, ring.add(t1.coefficient(m), t2.coefficient(m)));
		}
		return getPolynomial(coeff);
	}

	@Override
	public MultivariatePolynomial<T> negative(Polynomial<T> t) {
		Map<Monomial, T> coeff = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			coeff.put(m, this.ring.negative(t.coefficient(m)));
		}
		return getPolynomial(coeff);
	}

	@Override
	public MultivariatePolynomial<T> multiply(Polynomial<T> t1, Polynomial<T> t2) {
		SortedMap<Monomial, T> result = new TreeMap<>();
		for (Monomial m1 : t1.monomials()) {
			for (Monomial m2 : t2.monomials()) {
				Monomial m = m1.multiply(m2);
				T newValue = this.ring.multiply(t1.coefficient(m1), t2.coefficient(m2));
				if (result.containsKey(m)) {
					result.put(m, this.ring.add(result.get(m), newValue));
				} else {
					result.put(m, newValue);
				}
			}
		}
		return this.getPolynomial(result);
	}

	@Override
	public MultivariatePolynomial<T> multiply(T a, Polynomial<T> t2) {
		return this.multiply(a, getMonomial(new int[this.numvars]), t2);
	}

	public MultivariatePolynomial<T> multiply(T a, Monomial t1, Polynomial<T> t2) {
		Map<Monomial, T> coeff = new TreeMap<>();
		for (Monomial m : t2.monomials()) {
			coeff.put(m.multiply(t1), this.ring.multiply(a, t2.coefficient(m)));
		}
		return getPolynomial(coeff);
	}

	@Override
	public MultivariatePolynomial<T> one() {
		return this.getEmbedding(this.ring.one());
	}

	@Override
	public boolean isUnit(Polynomial<T> t) {
		return t.degree() == 0 && ring.isUnit(t.coefficient(getMonomial(new int[numvars])));
	}

	@Override
	public int numberOfVariables() {
		return this.numvars;
	}

	@Override
	public boolean isZeroDivisor(Polynomial<T> t) {
		throw new UnsupportedOperationException("Not implemented");
	}

	@Override
	public boolean isDivisible(Polynomial<T> dividend, Polynomial<T> divisor) {
		return remainder(dividend, divisor).equals(zero());
	}

	@Override
	public BigInteger euclidMeasure(Polynomial<T> t) {
		return BigInteger.valueOf(t.degree());
	}

	@Override
	public BigInteger getNumberOfUnits() {
		return this.ring.getNumberOfUnits();
	}

	@Override
	public Iterable<Polynomial<T>> getUnits() {
		return new Iterable<Polynomial<T>>() {
			private Iterable<T> elements = ring.getMultiplicativeGroup();

			@Override
			public Iterator<Polynomial<T>> iterator() {
				return new Iterator<Polynomial<T>>() {
					private Iterator<T> it = elements.iterator();

					@Override
					public boolean hasNext() {
						return it.hasNext();
					}

					@Override
					public MultivariatePolynomial<T> next() {
						T result = it.next();
						return MultivariatePolynomialRing.this.getEmbedding(result);
					}

				};
			}

		};
	}

	@Override
	public MultivariatePolynomial<T> inverse(Polynomial<T> t) {
		if (!this.isUnit(t))
			throw new RuntimeException();
		return this.getEmbedding(this.ring.inverse(t.coefficient(getMonomial(new int[this.numvars]))));
	}

	public Comparator<Monomial> getComparator() {
		return comparator;
	}

	@Override
	public Monomial getMonomial(int[] exponents) {
		return new Monomial(comparator, exponents);
	}

	@Override
	public T evaluate(Polynomial<T> t, List<T> ts) {
		MultivariatePolynomial<T> result = this.partiallyEvaluate(t, ts);
		if (result.degree() > 0)
			throw new ArithmeticException("not enough values");
		Monomial m = getMonomial(new int[ts.size()]);
		return result.coefficient(m);
	}

	@Override
	public MultivariatePolynomial<T> partiallyEvaluate(Polynomial<T> t, List<T> ts) {
		if (ts.size() != this.numvars)
			throw new ArithmeticException("wrong number of vars");
		Map<Monomial, T> coeff = new TreeMap<Monomial, T>();
		for (Monomial m : t.monomials()) {
			Monomial mp = this.dehomMonomialPartially(m, ts);
			if (!coeff.containsKey(mp))
				coeff.put(mp, this.ring.zero());
			coeff.put(mp, this.ring.add(coeff.get(mp),
					this.ring.multiply(t.coefficient(m), this.evaluateMonomialPartially(m, ts))));
		}
		return getPolynomial(coeff);
	}

	private T evaluateMonomialPartially(Monomial m, List<T> values) {
		T result = this.ring.one();
		for (int i = 0; i < this.numvars; i++)
			if (values.get(i) != null)
				result = this.ring.multiply(result, this.ring.power(values.get(i), m.exponents()[i]));
		return result;
	}

	private Monomial dehomMonomialPartially(Monomial m, List<T> values) {
		int[] exponents = new int[numvars];
		for (int i = 0; i < this.numvars; i++)
			if (values.get(i) == null)
				exponents[i] = m.exponents()[i];
		return getMonomial(exponents);
	}

	@Override
	public MultivariatePolynomial<T> substitute(Polynomial<T> t, List<Polynomial<T>> values) {
		MultivariatePolynomial<T> result = this.zero();
		for (Monomial m : t.monomials()) {
			MultivariatePolynomial<T> value = getEmbedding(t.coefficient(m));
			for (int i = 0; i < t.numberOfVariables(); i++) {
				value = multiply(value, power(values.get(i), m.exponents()[i]));
			}
			result = add(result, value);
		}
		return result;
	}

	@Override
	public boolean isHomogeneous(Polynomial<T> t) {
		for (Monomial m : t.monomials())
			if (m.degree() != t.degree())
				return false;
		return true;
	}

	@Override
	public boolean canBeNormalized(Polynomial<T> t) {
		T lc = t.leadingCoefficient();
		for (Monomial m : t.monomials()) {
			if (!ring.isDivisible(t.coefficient(m), lc)) {
				return false;
			}
		}
		return true;
	}

	public MultivariatePolynomial<T> normalize(Polynomial<T> t) {
		return divideScalar(t, t.leadingCoefficient());
	}

	@Override
	public T content(Polynomial<T> t) {
		T content = ring.zero();
		for (Monomial m : t.monomials()) {
			content = ring.gcd(t.coefficient(m), content);
		}
		return content;
	}

	@Override
	public MultivariatePolynomial<T> divideScalar(Polynomial<T> dividend, T divisor) {
		Map<Monomial, T> coeff = new TreeMap<Monomial, T>();
		for (Monomial m : dividend.monomials()) {
			coeff.put(m, ring.divide(dividend.coefficient(m), divisor));
		}
		return getPolynomial(coeff);
	}

	@Override
	public MultivariatePolynomial<T> contentFree(Polynomial<T> t) {
		return divideScalar(t, content(t));
	}

	@Override
	public MultivariatePolynomial<T> derivative(Polynomial<T> t, int variable) {
		if (variable < 1 || variable > this.numvars)
			throw new ArithmeticException("wrong number of vars");
		variable--;
		Map<Monomial, T> coeff = new TreeMap<Monomial, T>();
		for (Monomial m : t.monomials()) {
			int k = m.exponents()[variable];
			if (k == 0)
				continue;
			int[] exponents = new int[this.numvars];
			for (int i = 0; i < this.numvars; i++)
				exponents[i] = m.exponents()[i];
			exponents[variable]--;
			Monomial mp = getMonomial(exponents);
			if (!coeff.containsKey(mp))
				coeff.put(mp, ring.zero());
			coeff.put(mp, ring.add(coeff.get(mp), ring.multiply(k, t.coefficient(m))));
		}
		return getPolynomial(coeff);
	}

	@Override
	public Polynomial<T> resultant(Polynomial<T> t1, Polynomial<T> t2, int variable) {
		UnivariatePolynomialRing<Polynomial<T>> eliminatedRing = eliminateVariable().getUnivariatePolynomialRing();
		return eliminatedRing.resultant(asUnivariatePolynomial(t1, variable), asUnivariatePolynomial(t2, variable));
	}

	@Override
	public Polynomial<T> gcd(Polynomial<T> t1, Polynomial<T> t2) {
//		if (ring.equals(Rationals.q())) {
//			return gcdOverRationals(t1, t2);
//		}
		UnivariatePolynomialRing<Polynomial<T>> eliminatedRing = eliminateVariable().getUnivariatePolynomialRing();
		return fromUnivariatePolynomial(
				eliminatedRing.gcd(asUnivariatePolynomial(t1, 1), asUnivariatePolynomial(t2, 1)), 1);
	}

	private Polynomial<T> gcdOverRationals(Polynomial<T> t1, Polynomial<T> t2) {
		Rationals q = Rationals.q();
		Polynomial<Fraction> p1 = (Polynomial<Fraction>) t1;
		Polynomial<Fraction> p2 = (Polynomial<Fraction>) t2;
		return null;
	}

	@Override
	public MultivariatePolynomial<T> round(Polynomial<T> t, int degree) {
		Map<Monomial, T> map = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			if (m.degree() > degree) {
				continue;
			}
			map.put(m, t.coefficient(m));
		}
		return getPolynomial(map);
	}
}
