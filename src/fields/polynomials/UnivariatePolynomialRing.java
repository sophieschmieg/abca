package fields.polynomials;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.vectors.FreeModule;
import fields.vectors.Vector;

public class UnivariatePolynomialRing<T extends Element<T>> extends AbstractPolynomialRing<T>
		implements PolynomialRing<T> {
	private Ring<T> ring;
	private SortedSet<Monomial> monomials;
	private Map<UnivariatePolynomial<T>, List<Polynomial<T>>> squareFreeFactorizatonCache;

	public UnivariatePolynomialRing(Ring<T> ring) {
		super(false);
		this.ring = ring;
		this.monomials = new TreeSet<>();
		this.monomials.add(getMonomial(new int[] {0}));
		this.squareFreeFactorizatonCache = new TreeMap<>();
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
	public BigInteger characteristic() {
		return ring.characteristic();
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
	public boolean isIntegral() {
		return ring.isIntegral();
	}

	@Override
	public boolean isEuclidean() {
		return this.ring.isIntegral() && this.ring.krullDimension() == 0;
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public UnivariatePolynomial<T> zero() {
		return new UnivariatePolynomial<T>(this, Collections.emptyList());
	}

	@Override
	public UnivariatePolynomial<T> one() {
		return new UnivariatePolynomial<T>(this, Collections.singletonList(ring.one()));
	}

	public UnivariatePolynomial<T> getVar() {
		return getVarPower(1);
	}

	@Override
	public Polynomial<T> getVar(int i) {
		return getVar();
	}

	public UnivariatePolynomial<T> getVarPower(int power) {
		return getEmbedding(ring.one(), power);
	}

	@Override
	public UnivariatePolynomial<T> getVarPower(int i, int power) {
		return getVarPower(power);
	}

	@Override
	public UnivariatePolynomial<T> getEmbedding(T t) {
		return new UnivariatePolynomial<T>(this, Collections.singletonList(t));
	}

	public UnivariatePolynomial<T> getPolynomial(List<T> coefficients) {
		return new UnivariatePolynomial<T>(this, coefficients);
	}

	@SafeVarargs
	public final UnivariatePolynomial<T> getPolynomial(T... coefficients) {
		return getPolynomial(Arrays.asList(coefficients));
	}

	@Override
	public UnivariatePolynomial<T> getPolynomial(Map<Monomial, T> t) {
		List<T> coefficients = new ArrayList<>();
		for (Monomial m : t.keySet()) {
			int degree = m.exponents()[0];
			if (coefficients.size() <= degree) {
				for (int i = coefficients.size(); i <= degree; i++) {
					coefficients.add(ring.zero());
				}
			}
			coefficients.set(degree, t.get(m));
		}
		return new UnivariatePolynomial<T>(this, coefficients);
	}

	public UnivariatePolynomial<T> getEmbedding(T t, int exponent) {
		List<T> c = new ArrayList<>();
		for (int i = 0; i < exponent; i++) {
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

	public UnivariatePolynomial<T> toUnivariate(Polynomial<T> t) {
		if (t instanceof UnivariatePolynomial<?>) {
			return (UnivariatePolynomial<T>) t;
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
	public UnivariatePolynomial<T> getEmbeddingWithElimination(Polynomial<T> t, int shift) {
		if (t.numberOfVariables() + shift != 1) {
			throw new ArithmeticException("Number of variables don't match!");
		}
		List<T> values = new ArrayList<>();
		for (int i = 0; i < -shift; i++) {
			values.add(ring.one());
		}
		values.add(null);
		return toUnivariate(t.getPolynomialRing().partiallyEvaluate(t, values));
	}

	public UnivariatePolynomial<T> getEmbedding(Polynomial<T> t, int[] map) {
		Map<Monomial, T> c = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			int[] newexp = new int[1];
			for (int i = 0; i < t.numberOfVariables(); i++) {
				if (m.exponents()[i] != 0)
					newexp[map[i]] = m.exponents()[i];
			}
			c.put(getMonomial(newexp), t.coefficient(m));
		}
		return getPolynomial(c);
	}

	@Override
	public Polynomial<T> getRandomElement() {
		List<T> c = new ArrayList<>();
		int degree = new Random().nextInt(4);
		for (int i = 0; i <= degree; i++) {
			c.add(ring.getRandomElement());
		}
		return getPolynomial(c);
	}

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
	public UnivariatePolynomial<T> add(Polynomial<T> s1, Polynomial<T> s2) {
		UnivariatePolynomial<T> p1 = (UnivariatePolynomial<T>) s1;
		UnivariatePolynomial<T> p2 = (UnivariatePolynomial<T>) s2;
		List<T> c = new ArrayList<T>();
		int degree = Math.max(s1.degree(), s2.degree());
		for (int i = 0; i <= degree; i++) {
			c.add(ring.add(p1.univariateCoefficient(i), p2.univariateCoefficient(i)));
		}
		return getPolynomial(c);
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
	public UnivariatePolynomial<T> multiply(Polynomial<T> t1, Polynomial<T> t2) {
		UnivariatePolynomial<T> p1 = (UnivariatePolynomial<T>) t1;
		UnivariatePolynomial<T> p2 = (UnivariatePolynomial<T>) t2;
		List<T> c = new ArrayList<T>();
		int degree = p1.degree() + p2.degree();
		for (int i = 0; i <= degree; i++) {
			T a = ring.zero();
			for (int j = 0; j <= i; j++) {
				a = ring.add(a, ring.multiply(p1.univariateCoefficient(j), p2.univariateCoefficient(i - j)));
			}
			c.add(a);
		}
		return getPolynomial(c);
	}

	@Override
	public UnivariatePolynomial<T> multiply(T t1, Polynomial<T> t2) {
		return multiplyPower(0, t1, t2);
	}

	public UnivariatePolynomial<T> multiplyPower(int power, T t1, Polynomial<T> t2) {
		UnivariatePolynomial<T> p = (UnivariatePolynomial<T>) t2;
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
	public UnivariatePolynomial<T> divide(Polynomial<T> dividend, T divisor) {
		UnivariatePolynomial<T> p = (UnivariatePolynomial<T>) dividend;
		List<T> c = new ArrayList<T>();
		int degree = p.degree();
		for (int i = 0; i <= degree; i++) {
			c.add(ring.quotientAndRemainder(p.univariateCoefficient(i), divisor).get(0));
		}
		return getPolynomial(c);
	}

	@Override
	public boolean isZeroDivisor(Polynomial<T> t) {
		UnivariatePolynomial<T> p = (UnivariatePolynomial<T>) t;
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
	public UnivariatePolynomial<T> inverse(Polynomial<T> t) {
		return getEmbedding(ring.inverse(((UnivariatePolynomial<T>) t).univariateCoefficient(0)));
	}

	@Override
	public boolean isDivisible(Polynomial<T> dividend, Polynomial<T> divisor) {
		return quotientAndRemainder(dividend, divisor).get(1).equals(zero());
	}

	@Override
	public List<Polynomial<T>> generalQuotientAndRemainder(Polynomial<T> dividend, List<Polynomial<T>> divisors) {
		List<Polynomial<T>> result = new ArrayList<>();
		UnivariatePolynomial<T> remainder = toUnivariate(dividend);
		for (Polynomial<T> p : divisors) {
			UnivariatePolynomial<T> quotient = zero();
			while (remainder.degree() >= p.degree()
					&& ring.isDivisible(remainder.leadingCoefficient(), p.leadingCoefficient())) {
				int exponent = remainder.degree() - p.degree();
				T factor = ring.divide(remainder.leadingCoefficient(), p.leadingCoefficient());
				int deg = remainder.degree();
				remainder = toUnivariate(subtract(remainder, multiplyPower(exponent, factor, p)));
				if (remainder.degree() == deg) {
					System.err.println(remainder);
					System.err.println(factor);
					System.err.println(p);
					System.err.println(multiplyPower(exponent, factor, p));
				}
				quotient = add(quotient, getEmbedding(factor, exponent));
			}
			result.add(quotient);
		}
		result.add(remainder);
		return result;
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

	public List<Polynomial<T>> squareFreeFactorization(Polynomial<T> t) {
		UnivariatePolynomial<T> f = (UnivariatePolynomial<T>) t;
		if (this.squareFreeFactorizatonCache.containsKey(f)) {
			return this.squareFreeFactorizatonCache.get(f);
		}
		if (!t.leadingCoefficient().equals(this.ring.one())) {
			throw new RuntimeException("Not normalized.");
		}
		List<Polynomial<T>> result = new ArrayList<>();
		UnivariatePolynomial<T> tPrime = derivative(t);
		// Over algebraic closure, this decomposes into the factors of t with one degree
		// less.
		UnivariatePolynomial<T> degreeMinusI = normalize(gcd(t, tPrime));
		// Over the algebraic closure, this decomposes into the linear factors of t.
		UnivariatePolynomial<T> linearFactors = (UnivariatePolynomial<T>) this.divide(t, degreeMinusI);
		int i = 0;
		// Finding all factors of degree i with gcd(i, char F) = 1.
		while (linearFactors.degree() != 0) {
			i++;
			// Over the algebraic closure, this is all the factors with multiplicity higher
			// than i.
			UnivariatePolynomial<T> linearFactorsHigherDegree = (UnivariatePolynomial<T>) normalize(this.gcd(linearFactors,
					degreeMinusI));
			// This are all the linear factors that have multiplicity i in t.
			UnivariatePolynomial<T> factorThisDegree = (UnivariatePolynomial<T>) this.divide(linearFactors,
					linearFactorsHigherDegree);
			for (int j = 0; j < i; j++) {
				result.add(factorThisDegree);
			}
			// Updating the linear factors.
			linearFactors = linearFactorsHigherDegree;
			// Removing multiplicity i factors.
			degreeMinusI = (UnivariatePolynomial<T>) this.divide(degreeMinusI, linearFactorsHigherDegree);
		}
		// Remaining factors are multiplicity kp
		if (!degreeMinusI.equals(this.one())) {
			if (this.characteristic().equals(BigInteger.ZERO)) {
				return result;
			}
			// Substituting x^1/p.
			UnivariatePolynomial<T> pthRoot = this.zero();
			for (int m = f.degree(); m >= 0; m--) {
				if (m % this.characteristic().intValueExact() != 0 && !f.univariateCoefficient(m).equals(ring.zero())) {
					throw new ArithmeticException("SquareFreeFactorization wrong");
				}
				// x^q == x => x^(q/p)^p = x^(p^(n-1))^p = x.
				BigInteger exp = this.ring.getNumberOfElements().divide(this.ring.characteristic());
				T coeff = degreeMinusI.univariateCoefficient(m);
				pthRoot = this.add(pthRoot, this.multiply(this.ring.power(coeff, exp),
						this.getVarPower(1, m / this.characteristic().intValueExact())));
			}
			List<Polynomial<T>> factors = this.squareFreeFactorization(pthRoot);
			for (Polynomial<T> factor : factors) {
				for (int j = 0; j < this.characteristic().intValueExact(); j++) {
					result.add(factor);
				}
			}
		}
		this.squareFreeFactorizatonCache.put(f, result);
		return result;
	}

	@Override
	public List<Polynomial<T>> getAlgebraGenerators() {
		return Collections.singletonList(getVar());
	}

	@Override
	public boolean isLinearIndependent(List<Polynomial<T>> s) {
		int degree = 0;
		for (Polynomial<T> p : s) {
			degree = Math.max(degree, p.degree());
		}
		FreeModule<T> free = new FreeModule<>(ring, degree + 1);
		List<Vector<T>> asVectors = new ArrayList<>();
		for (Polynomial<T> p : s) {
			asVectors.add(asVector((UnivariatePolynomial<T>) p, degree));
		}
		return free.isLinearIndependent(asVectors);
	}

	@Override
	public boolean isGeneratingModule(List<Polynomial<T>> s) {
		return false;
	}

	@Override
	public List<Polynomial<T>> getModuleGenerators() {
		throw new InfinityException();
	}

	@Override
	public Vector<T> asVector(Polynomial<T> s) {
		return asVector(s, s.degree());
	}

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
	public Iterator<Polynomial<T>> iterator() {
		throw new InfinityException();
	}

	@Override
	public final T evaluate(Polynomial<T> t, List<T> ts) {
		if (ts.size() != 1) {
			throw new ArithmeticException("Only a single variable needed");
		}
		UnivariatePolynomial<T> p = (UnivariatePolynomial<T>) t;
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
	public UnivariatePolynomial<T> substitute(Polynomial<T> t, List<Polynomial<T>> values) {
		UnivariatePolynomial<T> p = (UnivariatePolynomial<T>) t;
		UnivariatePolynomial<T> result = this.zero();
		for (int i = 0; i <= t.degree(); i++) {
			result = add(result, multiply(p.univariateCoefficient(i), power(values.get(0), i)));
		}
		return result;
	}

	@Override
	public boolean canBeNormalized(Polynomial<T> t) {
		return ring.isDivisible(content(t), t.leadingCoefficient());
	}

	@Override
	public UnivariatePolynomial<T> normalize(Polynomial<T> t) {
		return divide(t, t.leadingCoefficient());
	}

	@Override
	public T content(Polynomial<T> t) {
		T content = ring.zero();
		UnivariatePolynomial<T> p = (UnivariatePolynomial<T>) t;
		for (int i = 0; i <= p.degree(); i++) {
			content = ring.gcd(content, p.univariateCoefficient(i));
		}
		return content;
	}

	@Override
	public UnivariatePolynomial<T> contentFree(Polynomial<T> t) {
		return divide(t, content(t));
	}

	@Override
	public UnivariatePolynomial<T> derivative(Polynomial<T> t, int variable) {
		return derivative(t);
	}

	public UnivariatePolynomial<T> derivative(Polynomial<T> t) {
		UnivariatePolynomial<T> p = (UnivariatePolynomial<T>) t;
		List<T> c = new ArrayList<>();
		for (int i = 0; i < p.degree(); i++) {
			c.add(ring.multiply(i + 1, p.univariateCoefficient(i + 1)));
		}
		return getPolynomial(c);
	}

	public T trace(Polynomial<T> t, int degree) {
		return ((UnivariatePolynomial<T>) t).trace(degree);
	}
	
	@Override
	public UnivariatePolynomial<T> round(Polynomial<T> t, int degree) {
		UnivariatePolynomial<T> p = toUnivariate(t);
		return getPolynomial(p.coefficients().subList(0, Math.min(degree, p.degree()) + 1));
	}
}
