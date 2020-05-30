package fields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.CoordinateRing.CoordinateRingElement;
import fields.Polynomial.EliminationOrder;
import fields.Polynomial.Monomial;
import util.Pair;

public class PolynomialRing<T extends Element> extends AbstractAlgebra<T, Polynomial<T>> {
	private Field<T> field;
	private int numvars;
	private Comparator<Polynomial.Monomial> comparator;
	private Map<Polynomial<T>, List<Polynomial<T>>> squareFreeFactorizatonCache;
	private Map<Polynomial<T>, List<Pair<Polynomial<T>, Integer>>> distinctDegreeFactorizatonCache;
	private Map<Polynomial<T>, Integer> distinctDegreeFactorizationMaxDegree;
	private Map<Polynomial<T>, Map<Integer, List<Polynomial<T>>>> irreducibleFactorizatonCache;

	public PolynomialRing(Field<T> field, int numvars, Comparator<Polynomial.Monomial> comparator) {
		this.field = field;
		this.numvars = numvars;
		this.comparator = comparator;
		this.squareFreeFactorizatonCache = new TreeMap<>();
		this.distinctDegreeFactorizatonCache = new TreeMap<>();
		this.distinctDegreeFactorizationMaxDegree = new TreeMap<>();
		this.irreducibleFactorizatonCache = new TreeMap<>();		
	}
	public Field<T> getField() {
		return this.field;
	}
	
	public Ring<T> getRing() {
		return this.field;
	}
	
	public List<Polynomial<T>> getGenerators() {
		List<Polynomial<T>> result = new ArrayList<>();
		for (int i = 1; i <= numvars; i++) {
			result.add(getVar(i));
		}
		return result;
	}

	public Polynomial<T> getEmbedding(T t, int[] exponents) {
		return this.getEmbedding(t, new Polynomial.Monomial(exponents));
	}
	public Polynomial<T> getEmbedding(T t, Polynomial.Monomial monomial) {
		Map<Polynomial.Monomial,T> c = new TreeMap<Polynomial.Monomial, T>(comparator);
		c.put(monomial, t);
		return new Polynomial<T>(c, this);
	}
	public Polynomial<T> getEmbedding(T t) {
		return this.getEmbedding(t, new int[numvars]);
	}
	public Polynomial<T> getEmbedding(Polynomial<T> t, int[] map) {
		Map<Polynomial.Monomial,T> c = new TreeMap<Polynomial.Monomial, T>(comparator);
		for (Polynomial.Monomial m : t.getCoefficients().keySet()) {
			int[] newexp = new int[this.numvars];
			for (int i = 0; i < t.getNumVars(); i++) {
				if (m.getExponents()[i] != 0)
					newexp[map[i]] = m.getExponents()[i];
			}
			c.put(new Polynomial.Monomial(newexp), t.getCoefficient(m));
		}
		return new Polynomial<T>(c, this);
	}
	public Polynomial<T> getEmbeddingShift(Polynomial<T> t, int shift) {
		if (this.getNumVars() != t.getNumVars() + shift)
			throw new ArithmeticException("Number of variables don't match!");
		int[] map = new int[t.getNumVars()];
		for (int i = 0; i < t.getNumVars(); i++)
			map[i] = i + shift;
		return this.getEmbedding(t, map);
	}
	public Polynomial<T> getPolynomial(Map<Polynomial.Monomial,T> coeff) {
		return new Polynomial<T>(coeff, this);
	}
	public Polynomial<T> getVar(int var) {
		return this.getVar(var, 1);
	}
	public Polynomial<T> getVar(int var, int power) {
		int[] exponents = new int[numvars];
		if (var < 1 || var > numvars)
			throw new ArithmeticException("wrong number of vars");
		var--;
		exponents[var] = power;
		return this.getEmbedding(this.field.one(), exponents);
	}

	public Polynomial<T> getLinear(List<T> coeff) {
		if (this.numvars != coeff.size())
			throw new ArithmeticException("size mismatch!");
		Map<Monomial, T> coeffs = new TreeMap<Monomial, T>(comparator);
		int[] exponents = new int[numvars];
		for (int i = 0; i < numvars; i++) {
			exponents[i] = 1;
			if (i > 0)
				exponents[i-1] = 0;
			coeffs.put(new Monomial(exponents), coeff.get(i));
		}
		return this.getPolynomial(coeffs);
	}
	@SafeVarargs
	public final Polynomial<T> getLinear(T... coeffs) {
		return this.getLinear(Arrays.asList(coeffs));
	}

	public PolynomialRing<T> addVariableWithElimination(int shift) {
		if (shift < 0)
			throw new ArithmeticException("Cannot count");
		return new PolynomialRing<T>(field, this.numvars + shift, new EliminationOrder(Polynomial.LEX, this.comparator, shift));
	}

	@Override
	public Polynomial<T> getRandomElement() {
		int degree = (int)Math.abs(new Random().nextGaussian() * 2.0);
		return getRandomElement(degree);
	}
	
	private Polynomial<T> getRandomElement(int degree, int var) {
		if (degree < 0) {
			return this.zero();
		}
		if (var > numvars) {
			return this.getEmbedding(field.getRandomElement());
		}
		Polynomial<T> result = this.zero();
		for (int i = 0; i <= degree; i++) {
			result = this.add(result, this.multiply(this.getVar(var, i), this.getRandomElement(degree - i, var + 1)));
		}
		return result;
	}
	
	public Polynomial<T> getRandomElement(int degree) {
		return this.getRandomElement(degree, 1);
	}
	@Override
	public BigInteger getNumberOfElements() {
		return BigInteger.valueOf(-1);
	}
	@Override
	public Iterator<Polynomial<T>> iterator() {
		throw new InfinityException();
	}
	@Override
	public Polynomial<T> zero() {
		return this.getEmbedding(this.field.zero());
	}
	@Override
	public Polynomial<T> add(Polynomial<T> t1, Polynomial<T> t2) {
		Map<Polynomial.Monomial,T> coeff = new TreeMap<Polynomial.Monomial, T>(this.comparator);
		coeff.put(new Monomial(new int[this.numvars]), this.field.zero());
		Set<Polynomial.Monomial> monomeset = new TreeSet<Polynomial.Monomial>(this.comparator);
		monomeset.addAll(t1.getCoefficients().keySet());
		monomeset.addAll(t2.getCoefficients().keySet());
		for (Polynomial.Monomial m : monomeset) {
			T c = this.field.zero();
			if (t1.getCoefficients().containsKey(m))
				c = this.field.add(c, t1.getCoefficients().get(m));
			if (t2.getCoefficients().containsKey(m))
				c = this.field.add(c, t2.getCoefficients().get(m));
			coeff.put(m, c);
		}
		return new Polynomial<T>(coeff, this);
	}
	@Override
	public Polynomial<T> negative(Polynomial<T> t) {
		Map<Polynomial.Monomial,T> coeff = new TreeMap<Polynomial.Monomial, T>(this.comparator);
		coeff.put(new Monomial(new int[this.numvars]), this.field.zero());
		for (Polynomial.Monomial m : t.getCoefficients().keySet()) {
			coeff.put(m, this.field.negative(t.getCoefficients().get(m)));
		}
		return new Polynomial<T>(coeff, this);
	}
	@Override
	public Polynomial<T> multiply(Polynomial<T> t1, Polynomial<T> t2) {
		SortedMap<Monomial, T> result = new TreeMap<>(this.comparator);
		for (Polynomial.Monomial m1 : t1.getCoefficients().keySet()) {
			for (Polynomial.Monomial m2 : t2.getCoefficients().keySet()) {
				Monomial m = this.multiply(m1, m2);
				T newValue = this.field.multiply(t1.getCoefficient(m1), t2.getCoefficient(m2));
				if (result.containsKey(m)) {
					result.put(m, this.field.add(result.get(m), newValue));
				} else {
					result.put(m, newValue);
				}
			}
		}
		return this.getPolynomial(result);
	}
	
	public Polynomial<T> multiply(T a, Polynomial<T> t2) {
		return this.multiply(a, new Polynomial.Monomial(new int[this.numvars]), t2);
	}
	
	public Polynomial<T> scalarMultiply(T a, Polynomial<T> t) {
		return this.multiply(a, t);
	}
	
	public Polynomial<T> multiply(T a, Polynomial.Monomial t1, Polynomial<T> t2) {
		Map<Polynomial.Monomial, T> coeff = new TreeMap<Polynomial.Monomial, T>(this.comparator);
		coeff.put(new Monomial(new int[this.numvars]), this.field.zero());
		for (Polynomial.Monomial m : t2.getCoefficients().keySet()) {
			coeff.put(Polynomial.Monomial.multiply(m, t1), this.field.multiply(a, t2.getCoefficients().get(m)));
		}
		return new Polynomial<T>(coeff, this);
	}
	@Override
	public Polynomial<T> one() {
		return this.getEmbedding(this.field.one());
	}
	@Override
	public BigInteger characteristic() {
		return this.field.characteristic();
	}
	@Override
	public boolean isUnit(Polynomial<T> t) {
		return t.getDegree() == 0;
	}
	@Override
	public boolean isFinite() {
		return false;
	}
	
	public boolean isFree() {
		return true;
	}
	
	public int getNumVars() {
		return this.numvars;
	}
	public boolean hasRoots(Polynomial<T> t) {
		if (this.numvars != 1) {
			throw new RuntimeException("Not implemented.");
		}
		if (!this.field.isFinite()) {
			throw new RuntimeException("Not implemented.");
		}
		CoordinateRing<T> cr = new CoordinateRing<T>(this, this.getIdeal(Collections.singletonList(t)));
		CoordinateRingElement<T> x = cr.getEmbedding(this.getVar(1));
		CoordinateRingElement<T> xq = cr.power(x, this.field.getNumberOfElements());
		Polynomial<T> xqmx = cr.subtract(xq, x).getElement();
		List<Polynomial<T>> egcd = this.extendedEuclidean(t, xqmx);
		return egcd.get(0).getDegree() > 0;
	}
	public Polynomial<T> gcd(Polynomial<T> t1, Polynomial<T> t2) {
		return this.extendedEuclidean(t1, t2).get(0).normalize();
	}
	/**
	 * Square free factorization of an univariate polynomial over a finite field.
	 * @param t Any normalized univariate polynomial over a finite field.
	 * @return list of square free factors of t.
	 */
	// Output: List of square free factors of t.
	public List<Polynomial<T>> squareFreeFactorization(Polynomial<T> t) {
		if (this.squareFreeFactorizatonCache.containsKey(t)) {
			return this.squareFreeFactorizatonCache.get(t);
		}
		if (this.numvars != 1) {
			throw new RuntimeException("Not univariate.");
		}
		if (!this.field.isFinite()) {
			throw new RuntimeException("Field is not finite.");
		}
		if (!t.getLeadingCoefficient().equals(this.field.one())) {
			throw new RuntimeException("Not normalized.");
		}
		List<Polynomial<T>> result = new ArrayList<>();
		Polynomial<T> tPrime = t.derivative(1);
		// Over algebraic closure, this decomposes into the factors of t with one degree less.
		Polynomial<T> degreeMinusI = this.gcd(t, tPrime);
		// Over the algebraic closure, this decomposes into the linear factors of t.
		Polynomial<T> linearFactors = this.quotientAndRemainder(t, degreeMinusI).get(0);
		int i = 0;
		// Finding all factors of degree i with gcd(i, char F) = 1.
		while (linearFactors.getDegree() != 0) {
			i++;
			// Over the algebraic closure, this is all the factors with multiplicity higher than i.
			Polynomial<T> linearFactorsHigherDegree = this.gcd(linearFactors, degreeMinusI);
			// This are all the linear factors that have multiplicity i in t.
			Polynomial<T> factorThisDegree = this.quotientAndRemainder(linearFactors, linearFactorsHigherDegree).get(0);
			for (int j = 0; j < i; j++) {
				result.add(factorThisDegree);
			}
			// Updating the linear factors.
			linearFactors = linearFactorsHigherDegree;
			// Removing multiplicity i factors.
			degreeMinusI = this.quotientAndRemainder(degreeMinusI, linearFactorsHigherDegree).get(0);
		}
		// Remaining factors are multiplicity kp
		if (!degreeMinusI.equals(this.one())) {
			// Substituting x^1/p.
			Polynomial<T> pthRoot = this.zero();
			for (Monomial m : degreeMinusI.getCoefficients().keySet()) {
				if (m.degree() % this.characteristic().intValueExact() != 0) {
					throw new ArithmeticException("SquareFreeFactorization wrong");
				}
				// x^q == x => x^(q/p)^p = x^(p^(n-1))^p = x.
				BigInteger exp = this.field.getNumberOfElements().divide(this.field.characteristic());
				T coeff = degreeMinusI.getCoefficient(m);
				pthRoot = this.add(pthRoot, this.multiply(this.field.power(coeff, exp), this.getVar(1, m.degree() / this.characteristic().intValueExact())));
			}
			List<Polynomial<T>> factors = this.squareFreeFactorization(pthRoot);
			for (Polynomial<T> factor : factors) {
				for (int j = 0; j < this.characteristic().intValueExact(); j++) {
					result.add(factor);
				}
			}
		}
		this.squareFreeFactorizatonCache.put(t, result);
		return result;
	}
	/**
	 * Factorizes a square free polynomial into factors that of a distinct degree.
	 * @param t univariate, normalized, square free polynomial over a finite field.
	 * @param maxDegree maximum degree of results needed. Set to -1 for all factors.
	 * @return list of pairs of polynomials and integers, such that the polynomials 
	 * are factors of t, that each consists of irreducible factors of the the given degree.
	 */
	public List<Pair<Polynomial<T>, Integer>> distinctDegreeFactorization(Polynomial<T> t, int maxDegree) {		
		if (this.distinctDegreeFactorizatonCache.containsKey(t)) {
			int previousMaxDegree = this.distinctDegreeFactorizationMaxDegree.get(t);
			if (previousMaxDegree == -1 || previousMaxDegree >= maxDegree) {
				return this.distinctDegreeFactorizatonCache.get(t);
			}
		}
		if (this.numvars != 1) {
			throw new RuntimeException("Not univariate.");
		}
		if (!this.field.isFinite()) {
			throw new RuntimeException("Field is not finite.");
		}
		if (!t.getLeadingCoefficient().equals(this.field.one())) {
			throw new RuntimeException("Not normalized.");
		}
		List<Pair<Polynomial<T>, Integer>> result = new ArrayList<>();		
		Polynomial<T> work = t;
		CoordinateRing<T> cr = new CoordinateRing<T>(this, this.getIdeal(Collections.singletonList(work)));
		CoordinateRingElement<T> x = cr.getEmbedding(this.getVar(1));
		CoordinateRingElement<T> xqi = x;		
		for (int degree = 1; work.getDegree() >= 2 * degree; degree++) {
			if (maxDegree >= 0 && degree > maxDegree) {
				break;
			}
			xqi = cr.power(xqi, this.field.getNumberOfElements());			
			CoordinateRingElement<T> xqimx = cr.subtract(xqi, x);
			Polynomial<T> gcd = this.gcd(work, xqimx.getElement());
			if (!gcd.equals(this.one())) {
				result.add(new Pair<>(gcd, degree));
				work = this.quotientAndRemainder(work, gcd).get(0);
				cr = new CoordinateRing<T>(this, this.getIdeal(Collections.singletonList(work)));
				x = cr.getEmbedding(this.getVar(1));
				xqi = x;
				for (int j = 0; j < degree; j++) {
					xqi = cr.power(xqi, this.field.getNumberOfElements());
				}
			}		
		}
		if (!work.equals(this.one())) {
			result.add(new Pair<>(work, work.getDegree()));
		}
		this.distinctDegreeFactorizatonCache.put(t, result);
		this.distinctDegreeFactorizationMaxDegree.put(t, maxDegree);
		return result;
	}
	/**
	 * Computes the irreducible factors of t, provided that t is univariate, normalized, square free,
	 * and only has irreducible factors of the given degree.
	 * @param t the polynomial to factorize.
	 * @param degree the degree of all irreducible factors.
	 * @return a list of irreducible factors of t.
	 */
	public List<Polynomial<T>> irreducibleFactorization(Polynomial<T> t, int degree) {
		if (!this.irreducibleFactorizatonCache.containsKey(t)) {
			this.irreducibleFactorizatonCache.put(t, new TreeMap<>());
		}
		Map<Integer, List<Polynomial<T>>> cache = this.irreducibleFactorizatonCache.get(t);
		if (cache.containsKey(degree)) {
			return cache.get(degree);
		}
		if (this.numvars != 1) {
			throw new RuntimeException("Not univariate.");
		}
		if (!this.field.isFinite()) {
			throw new RuntimeException("Field is not finite.");
		}
		if (!t.getLeadingCoefficient().equals(this.field.one())) {
			throw new RuntimeException("Not normalized.");
		}
		if (t.getDegree() % degree != 0) {
			throw new RuntimeException("Inputs impossible!");
		}
		int numberOfFactors = t.getDegree() / degree;
		CoordinateRing<T> cr = new CoordinateRing<T>(this, this.getIdeal(Collections.singletonList(t)));
		List<Polynomial<T>> factors = new ArrayList<>();
		factors.add(t);
		while (factors.size() < numberOfFactors) {
			// Generate random polynomial with degree smaller than t.
			Polynomial<T> h = this.zero();
			for (int i = 0; i < t.getDegree(); i++) {
				h = this.add(h, this.multiply(this.field.getRandomElement(), this.getVar(1, i)));
			}
			BigInteger one = BigInteger.ONE;
			BigInteger exponent = this.field.getNumberOfElements().pow(degree).subtract(one).shiftRight(1);
			Polynomial<T> g = cr.subtract(cr.power(cr.getEmbedding(h), exponent), cr.one()).getElement();
			List<Polynomial<T>> newFactors = new ArrayList<Polynomial<T>>();
			for (Polynomial<T> factor : factors) {
				if (factor.getDegree() == degree) {
					newFactors.add(factor);
					continue;
				}
				Polynomial<T> gcd = this.gcd(factor, g);
				if (!gcd.equals(this.one()) && !gcd.equals(factor)) {
					newFactors.add(gcd);
					newFactors.add(this.quotientAndRemainder(factor, gcd).get(0));
				} else {
					newFactors.add(factor);
				}
			}
			factors = newFactors;
		}
		cache.put(degree, factors);
		return factors;
	}
	public List<Polynomial<T>> factorization(Polynomial<T> t) {
		if (this.numvars != 1) {
			throw new RuntimeException("Multivariate factorization not implemented.");
		}
		if (!this.field.isFinite()) {
			throw new RuntimeException("Characteristic 0 factorization not implemented.");
		}
		List<Polynomial<T>> result = new ArrayList<>();
		if (!t.getLeadingCoefficient().equals(this.field.one())) {
			result.add(this.getEmbedding(t.getLeadingCoefficient()));
			t = t.normalize();
		}
		List<Polynomial<T>> squareFreeFactors = this.squareFreeFactorization(t);
		for (Polynomial<T> sff : squareFreeFactors) {
			List<Pair<Polynomial<T>, Integer>> distinctDegreeFactors = this.distinctDegreeFactorization(sff, -1);
			for (Pair<Polynomial<T>, Integer> ddf : distinctDegreeFactors) {
				result.addAll(this.irreducibleFactorization(ddf.getFirst(), ddf.getSecond()));
			}
		}
		return result;
	}
	public List<T> roots(Polynomial<T> t) {
		if (this.numvars != 1) {
			throw new RuntimeException("Multivariate root finding not implemented. Use ideal logic.");
		}
		if (!this.field.isFinite()) {
			throw new RuntimeException("Characteristic 0 root finding not implemented.");
		}
		List<T> result = new ArrayList<T>();
		t = t.normalize();
		List<Polynomial<T>> squareFreeFactors = this.squareFreeFactorization(t);
		for (Polynomial<T> sff : squareFreeFactors) {
			List<Pair<Polynomial<T>, Integer>> distinctDegreeFactors = this.distinctDegreeFactorization(sff, 1);
			for (Pair<Polynomial<T>, Integer> ddf : distinctDegreeFactors) {
				if (ddf.getSecond() != 1) {
					continue;
				}
				List<Polynomial<T>> linear = this.irreducibleFactorization(ddf.getFirst(), ddf.getSecond());
				for (Polynomial<T> linearFactor : linear) {
					result.add(this.field.negative(linearFactor.getCoefficient(new Monomial(new int[] {0}))));
				}
			}
		}
		return result;
	}
	@Override
	public boolean isIntegral() {
		return true;
	}
	@Override
	public boolean isZeroDivisor(Polynomial<T> t) {
		return t.equals(zero());
	}
	@Override
	public boolean isEuclidean() {
		return this.numvars == 1;
	}
	@Override
	public boolean isDivisible(Polynomial<T> dividend, Polynomial<T> divisor) {
		return this.quotientAndRemainder(dividend, divisor).get(1).equals(zero());
	}
	@Override
	public List<Polynomial<T>> quotientAndRemainder(Polynomial<T> dividend, Polynomial<T> divisor) {
		if (divisor.getDegree() > dividend.getDegree()) {
			List<Polynomial<T>> result = new ArrayList<>();
			result.add(this.zero());
			result.add(dividend);
			return result;
		}
		return this.generalQuotientAndRemainder(dividend, Collections.singletonList(divisor));
	}
	@Override
	public BigInteger euclidMeasure(Polynomial<T> t) {
		return BigInteger.valueOf(t.getDegree());
	}
	@Override
	public BigInteger getNumberOfUnits() {
		return this.field.getNumberOfElements().subtract(BigInteger.ONE);
	}
	@Override
	public Iterable<Polynomial<T>> getUnits() {
		return new Iterable<Polynomial<T>>() {
			private Iterable<T> elements = field.getMultiplicativeGroup();
			@Override
			public Iterator<Polynomial<T>> iterator() {
				return new Iterator<Polynomial<T>>() {
					private Iterator<T> it = elements.iterator();

					@Override
					public boolean hasNext() {
						return it.hasNext();
					}

					@Override
					public Polynomial<T> next() {
						T result = it.next();
						return PolynomialRing.this.getEmbedding(result);
					}

					@Override
					public void remove() {
						it.remove();
					}
				};
			}

		};
	}
	@Override
	public Polynomial<T> inverse(Polynomial<T> t) {
		if (!this.isUnit(t))
			throw new RuntimeException();
		return this.getEmbedding(this.field.inverse(t.getCoefficient(new Monomial(new int[this.numvars]))));
	}
	public Comparator<Polynomial.Monomial> getComparator() {
		return comparator;
	}
	public Ideal getIdeal(List<Polynomial<T>> generators) {
		return new Ideal(generators);
	}
	public SortedSet<Polynomial<T>> buchberger(List<Polynomial<T>> generators) {
		List<Polynomial<T>> groebner = new ArrayList<Polynomial<T>>();
		for (Polynomial<T> basispolynomial : generators) {
			Polynomial<T> reduced = this.reduce(basispolynomial, groebner);
			if (!reduced.equals(this.zero()))
				groebner.add(reduced.normalize());
		}
		for (int i = 1; i < groebner.size(); i++) {
			for (int j = 0; j < i; j++) { 
				Polynomial<T> pi = groebner.get(i);
				Polynomial<T> pj = groebner.get(j);
				Polynomial.Monomial li = pi.getLeadingTerm();
				Polynomial.Monomial lj = pj.getLeadingTerm();
				Polynomial.Monomial lcm = this.lcm(li, lj);
				li = this.divide(lcm, li);
				lj = this.divide(lcm, lj);
				Polynomial<T> pli = this.getEmbedding(this.field.one(), li); 
				Polynomial<T> plj = this.getEmbedding(this.field.one(), lj); 
				Polynomial<T> newGenerator = this.subtract(this.multiply(pli, pi), this.multiply(plj, pj));
				newGenerator = this.reduce(newGenerator, groebner);
				if (!newGenerator.equals(this.zero()))
					groebner.add(newGenerator.normalize());
			}
		}
		return this.reduceBasis(groebner);
	}
	public SortedSet<Polynomial<T>> reduceBasis(List<Polynomial<T>> generators) {
		List<Polynomial<T>> list = new ArrayList<Polynomial<T>>();
		list.addAll(generators);
		SortedSet<Polynomial<T>> reduced = new TreeSet<Polynomial<T>>(Collections.reverseOrder());
		for (int i = 0; i < generators.size(); i++) {
			list.remove(i);
			Polynomial<T> reduce = this.reduce(generators.get(i), list);
			list.add(i, generators.get(i));
			if (!reduce.equals(this.zero()))
				reduced.add(reduce.normalize());
		}
		return reduced;
	}
	public Polynomial<T> reduce(Polynomial<T> polynomial, Collection<Polynomial<T>> basis) {
		if (basis.isEmpty()) {
			return polynomial;
		}
		return this.generalQuotientAndRemainder(polynomial, new ArrayList<Polynomial<T>>(basis)).get(basis.size());
	}
	public List<Polynomial<T>> generalQuotientAndRemainder(Polynomial<T> polynomial, List<Polynomial<T>> basis) {
		List<SortedMap<Monomial, T>> quotient = new ArrayList<>();
		SortedMap<Monomial, T> remainder = new TreeMap<>(this.comparator);
		SortedMap<Monomial, T> p = new TreeMap<>(this.comparator);
		p.putAll(polynomial.getCoefficients());
		for (int i = 0; i < basis.size(); i++)
			quotient.add(new TreeMap<>(this.comparator));
		while (true) {
			if (p.isEmpty()) {
				List<Polynomial<T>> list = new ArrayList<Polynomial<T>>();
				for (SortedMap<Monomial, T> q : quotient) {
					list.add(this.getPolynomial(q));
				}
				list.add(this.getPolynomial(remainder));
				return list;
			}
			Polynomial.Monomial leadingmonomial = p.lastKey();
			T leadingcoefficient = p.get(leadingmonomial);
			boolean success = false;
			boolean totalFailure = true;
			for (int i = 0; i < basis.size(); i++) {
				Polynomial<T> base = basis.get(i);
				if (this.comparator.compare(leadingmonomial, base.getLeadingTerm()) >= 0) {
					totalFailure = false;
				}
				Polynomial.Monomial div = this.divide(leadingmonomial, base.getLeadingTerm());
				if (div == null)
					continue;
				T factor = this.field.divide(leadingcoefficient, base.getLeadingCoefficient());
				quotient.get(i).put(div, factor);
				for (Monomial m : base.getCoefficients().keySet()) {
					T subtrahend = this.field.multiply(factor, base.getCoefficient(m));
					Monomial multiplied = this.multiply(m, div);
					T newValue;
					if (p.containsKey(multiplied)) {
						newValue = this.field.subtract(p.get(multiplied), subtrahend);
					} else {
						newValue = this.field.negative(subtrahend);
					}
					if (newValue.equals(this.field.zero())) {
						p.remove(multiplied);
					} else {
						p.put(multiplied, newValue);
					}
				}
				success = true;
				break;
			}
			if (totalFailure) {
				remainder.putAll(p);
				p.clear();
			} else if (!success) {
				remainder.put(leadingmonomial, leadingcoefficient);
				p.remove(leadingmonomial);
			}
		}
	}
	private Polynomial.Monomial divide(Polynomial.Monomial m1, Polynomial.Monomial m2) {
		int[] dividedexponents = new int[this.numvars];
		for (int i = 0; i < this.numvars; i++) {
			if (m1.getExponents()[i] < m2.getExponents()[i])
				return null;
			dividedexponents[i] = m1.getExponents()[i] - m2.getExponents()[i];
		}
		return new Polynomial.Monomial(dividedexponents);
	}
	private Polynomial.Monomial multiply(Polynomial.Monomial m1, Polynomial.Monomial m2) {
		int[] exponents = new int[this.numvars];
		for (int i = 0; i < this.numvars; i++) {
			exponents[i] = m1.getExponents()[i] + m2.getExponents()[i];
		}
		return new Polynomial.Monomial(exponents);
	}
	private Polynomial.Monomial lcm(Polynomial.Monomial m1, Polynomial.Monomial m2) {
		int[] exponents = new int[this.numvars];
		for (int i = 0; i < this.numvars; i++) {
			exponents[i] = Math.max(m1.getExponents()[i], m2.getExponents()[i]);
		}
		return new Polynomial.Monomial(exponents);
	}
	public class Ideal extends AbstractModule<Polynomial<T>, Polynomial<T>> {
		private SortedSet<Polynomial<T>> basis;
		private int[][] leadingMonomials;
		private int[] set;
		private int dimension;

		public Ideal(List<Polynomial<T>> generators) {
			this.basis = PolynomialRing.this.buchberger(generators);
			this.leadingMonomials = new int[this.basis.size()][PolynomialRing.this.numvars];
			int i = 0;
			for (Polynomial<T> generator : this.basis) {
				this.leadingMonomials[i] = generator.getLeadingTerm().getExponents();
				i++;
			}
			this.set = new int[PolynomialRing.this.numvars];
			int[] tmp = this.makeset(0, new int[PolynomialRing.this.numvars + 1]);
			this.set = Arrays.copyOf(tmp, this.set.length);
			this.dimension = PolynomialRing.this.numvars - tmp[this.set.length];
		}
		public SortedSet<Polynomial<T>> getBasis() {
			return Collections.unmodifiableSortedSet(this.basis);
		}
		public boolean equalsIdeal(Ideal other) {
			return this.basis.equals(other.basis);
		}
		public boolean contains(Polynomial<T> polynomial) {
			return PolynomialRing.this.reduce(polynomial, new ArrayList<Polynomial<T>>(this.basis)).equals(PolynomialRing.this.zero());
		}
		public boolean contains(Ideal ideal) {
			for (Polynomial<T> gen : ideal.basis)
				if (!this.contains(gen))
					return false;
			return true;
		}
		public Polynomial<T> residue(Polynomial<T> f) {
			return PolynomialRing.this.reduce(f, this.basis);
		}
		public int dimension() {
			return this.dimension;
		}
		private int[] makeset(int level, int[] set) {
			if (level == this.basis.size())
				return set;
			int[] optimalset = new int[set.length];
			optimalset[set.length - 1] = -1;
			for (int i = 0; i < this.leadingMonomials[level].length; i++) {
				if (this.leadingMonomials[level][i] == 0)
					continue;
				if (set[i] == 1)
					return this.makeset(level + 1, set);
				//int[] myset = Arrays.copyOf(set, set.length);
				set[i] = 1;
				set[set.length - 1]++;
				int[] sethere = this.makeset(level + 1, set);
				if (optimalset[set.length - 1] == -1 || optimalset[set.length - 1] > sethere[set.length - 1])
					optimalset = Arrays.copyOf(sethere, sethere.length);
				set[i] = 0;
				set[set.length - 1]--;
			}
			return optimalset;
		}
		public Ideal add(Ideal ideal) {
			List<Polynomial<T>> generators = new ArrayList<Polynomial<T>>();
			generators.addAll(this.basis);
			generators.addAll(ideal.basis);
			return new Ideal(generators);
		}
		public Ideal multiply(Ideal ideal) {
			List<Polynomial<T>> generators = new ArrayList<Polynomial<T>>();
			for (Polynomial<T> p : this.basis) {
				for (Polynomial<T> q : ideal.basis) {
					generators.add(PolynomialRing.this.multiply(p, q));
				}
			}
			return new Ideal(generators);
		}
		public Ideal intersect(Ideal ideal) {
			List<Polynomial<T>> generators = new ArrayList<Polynomial<T>>();
			PolynomialRing<T> ring = PolynomialRing.this.addVariableWithElimination(1);
			Polynomial<T> t = ring.getVar(1);
			Polynomial<T> omt = ring.subtract(ring.one(), t);
			for (Polynomial<T> f : this.basis) {
				generators.add(ring.multiply(t, ring.getEmbeddingShift(f, 1)));
			}
			for (Polynomial<T> g : ideal.basis) {
				generators.add(ring.multiply(omt, ring.getEmbeddingShift(g, 1)));
			}
			Ideal intersectionIdeal = ring.getIdeal(generators);
			List<Polynomial<T>> intersectionGenerators = new ArrayList<Polynomial<T>>();
			for (Polynomial<T> b : intersectionIdeal.getBasis()) {
				if (b.getLeadingTerm().getExponents()[0] == 0)
					intersectionGenerators.add(PolynomialRing.this.getEmbeddingShift(b, -1));
			}
			return new Ideal(intersectionGenerators);
		}
		public Ideal saturate(Polynomial<T> by) {
			PolynomialRing<T> ring = PolynomialRing.this.addVariableWithElimination(1);
			Polynomial<T> t = ring.getVar(1);
			Polynomial<T> f = ring.getEmbeddingShift(by, 1);
			Polynomial<T> inv = ring.subtract(ring.one(), ring.multiply(f, t));
			List<Polynomial<T>> generators = new ArrayList<Polynomial<T>>();
			generators.add(inv);
			for (Polynomial<T> p : this.basis)
				generators.add(ring.getEmbeddingShift(p, 1));
			Ideal saturationIdeal = ring.getIdeal(generators);
			List<Polynomial<T>> saturationGenerators = new ArrayList<Polynomial<T>>();
			for (Polynomial<T> b : saturationIdeal.getBasis()) {
				if (b.getLeadingTerm().getExponents()[0] == 0)
					saturationGenerators.add(PolynomialRing.this.getEmbeddingShift(b, -1));
			}
			return new Ideal(saturationGenerators);

		}
		
		public String toString() {
			return this.basis.toString();
		}
		
		public Ring<Polynomial<T>> getRing() {
			return PolynomialRing.this;
		}
		
		public boolean isFree() {
			return true;
		}
		
		public boolean isFinite() {
			return false;
		}
		
		public Polynomial<T> zero() {
			return PolynomialRing.this.zero();
		}
		
		public Polynomial<T> add(Polynomial<T> t1, Polynomial<T> t2) {
			return PolynomialRing.this.add(t1, t2);
		}
		
		public Polynomial<T> negative(Polynomial<T> t) {
			return PolynomialRing.this.negative(t);
		}
		
		public Polynomial<T> scalarMultiply(Polynomial<T> t1, Polynomial<T> t2) {
			return PolynomialRing.this.multiply(t1, t2);
		}
		
		@Override
		public Polynomial<T> getRandomElement() {
			Polynomial<T> randomElement = this.zero();
			for (Polynomial<T> b : this.basis) {
				randomElement = this.add(randomElement, this.scalarMultiply(PolynomialRing.this.getRandomElement(), b));
			}
			return randomElement;
		}
		
		@Override
		public BigInteger getNumberOfElements() {
			return BigInteger.valueOf(-1);
		}
		@Override
		public Iterator<Polynomial<T>> iterator() {
			throw new InfinityException();
		}
	}
}
