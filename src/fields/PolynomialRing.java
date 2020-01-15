package fields;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.Polynomial.EliminationOrder;
import fields.Polynomial.Monomial;

public class PolynomialRing<T extends Element> extends AbstractRing<Polynomial<T>> {
	private Field<T> field;
	private int numvars;
	private Comparator<Polynomial.Monomial> comparator;
	
	public PolynomialRing(Field<T> field, int numvars, Comparator<Polynomial.Monomial> comparator) {
		this.field = field;
		this.numvars = numvars;
		this.comparator = comparator;
	}
	public Field<T> getField() {
		return this.field;
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
	/*public Ideal getLine(ProjectivePoint<T> p, ProjectivePoint<T> q) {
		if (p.equals(q))
			throw new ArithmeticException("Need two points for a line");
		Field<T> rh = field;
		List<T> list = new ArrayList<T>();
		list.add(rh.subtract(rh.multiply(p.getCoord(2), q.getCoord(3)), rh.multiply(p.getCoord(3), q.getCoord(2))));
		list.add(rh.subtract(rh.multiply(p.getCoord(3), q.getCoord(1)), rh.multiply(p.getCoord(1), q.getCoord(3))));
		list.add(rh.subtract(rh.multiply(p.getCoord(1), q.getCoord(2)), rh.multiply(p.getCoord(2), q.getCoord(1))));
		return Polynomial.getLine(rh, list, comparator);
	}
*/

	public PolynomialRing<T> addVariableWithElimination(int shift) {
		if (shift < 0)
			throw new ArithmeticException("Cannot count");
		return new PolynomialRing<T>(field, this.numvars + shift, new EliminationOrder(Polynomial.LEX, this.comparator, shift));
	}

	@Override
	public Polynomial<T> getRandomElement() {
		throw new UnsupportedOperationException();
	}
	@Override
	public int getNumberOfElements() {
		return -1;
	}
	@Override
	public Iterable<Polynomial<T>> getElements() throws InfinityException {
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
		Polynomial<T> result = this.zero();
		for (Polynomial.Monomial m : t1.getCoefficients().keySet()) 
			result = this.add(result, this.multiply(t1.getCoefficients().get(m), m, t2));
		return result;
	}
	public Polynomial<T> multiply(T a, Polynomial<T> t2) {
		return this.multiply(a, new Polynomial.Monomial(new int[this.numvars]), t2);
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
	public int characteristic() {
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
	public int getNumVars() {
		return this.numvars;
	}
	@Override
	public boolean isIntegral() {
		return true;
	}
	@Override
	public boolean isEuclidean() {
		return this.numvars == 1;
	}
	@Override
	public List<Polynomial<T>> quotientAndRemainder(Polynomial<T> dividend, Polynomial<T> divisor) {
		return this.generalQuotientAndRemainder(dividend, Collections.singletonList(divisor));
		/*
		Polynomial<T> quotient, remainder;
		List<Polynomial<T>> list;
		if (dividend.getDegree() < divisor.getDegree()) {
			list = new ArrayList<Polynomial<T>>();
			list.add(this.zero());
			list.add(dividend);
			return list;
		}
		T highdividend = this.getCoefficient(dividend, dividend.getDegree()); 
		T highdivisor = this.getCoefficient(divisor, divisor.getDegree());
		List<T> qr = this.ring.quotientAndRemainder(highdividend, highdivisor);
		if (!qr.get(1).equals(this.ring.zero()))
			throw new ArithmeticException("Not divisible");
		T factor = qr.get(0);
		Polynomial<T> pfactor = Polynomial.getEmbedding(this.ring, factor, 1, this.comparator);
		pfactor = this.multiply(pfactor, Polynomial.getVar(this.ring, 1, 1, dividend.getDegree() - divisor.getDegree(), this.comparator));
		dividend = this.add(dividend, this.multiply(this.negative(pfactor), divisor));
		list = this.quotientAndRemainder(dividend, divisor);
		quotient = this.add(pfactor, list.get(0));
		remainder = list.get(1);
		list = new ArrayList<Polynomial<T>>();
		list.add(quotient);
		list.add(remainder);
		return list;
	}
	private T getCoefficient(Polynomial<T> poly, int i) {
		int[] exponents = new int[] {i};
		return poly.getCoefficients().get(new Monomial(exponents));*/
	}
	@Override
	public int getNumberOfUnits() {
		return 0;
	}
	@Override
	public Iterable<Polynomial<T>> getUnits() {
		return new Iterable<Polynomial<T>>() {
			private Iterable<T> elements = field.getMultiplicativeGroup().getElements();
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
		return this.generalQuotientAndRemainder(polynomial, new ArrayList<Polynomial<T>>(basis)).get(basis.size());
	}
	public List<Polynomial<T>> generalQuotientAndRemainder(Polynomial<T> polynomial, List<Polynomial<T>> basis) {
		List<Polynomial<T>> quotient = new ArrayList<Polynomial<T>>();
		Polynomial<T> remainder = this.zero();
		for (int i = 0; i < basis.size(); i++)
			quotient.add(this.zero());
		while (true) {
			if (polynomial.equals(this.zero())) {
				List<Polynomial<T>> list = new ArrayList<Polynomial<T>>();
				list.addAll(quotient);
				list.add(remainder);
				return list;
			}
			Polynomial.Monomial leadingmonomial = polynomial.getLeadingTerm();
			T leadingcoefficient = polynomial.getLeadingCoefficient();
			boolean success = false;
			for (int i = 0; i < basis.size(); i++) {
				Polynomial.Monomial div = this.divide(leadingmonomial, basis.get(i).getLeadingTerm());
				if (div == null)
					continue;
				Polynomial<T> divpoly = this.getEmbedding(this.field.divide(leadingcoefficient, basis.get(i).getLeadingCoefficient()), div);
				quotient.set(i, this.add(divpoly, quotient.get(i)));
				polynomial = this.subtract(polynomial, this.multiply(divpoly, basis.get(i)));
				success = true;
				break;
			}
			if (!success) {
				Polynomial<T> leading = this.getEmbedding(leadingcoefficient, leadingmonomial); 
				remainder = this.add(remainder, leading);
				polynomial = this.subtract(polynomial, leading);
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
	private Polynomial.Monomial lcm(Polynomial.Monomial m1, Polynomial.Monomial m2) {
		int[] exponents = new int[this.numvars];
		for (int i = 0; i < this.numvars; i++) {
			exponents[i] = Math.max(m1.getExponents()[i], m2.getExponents()[i]);
		}
		return new Polynomial.Monomial(exponents);
	}
	public class Ideal {
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
			PolynomialRing<T> ring = PolynomialRing.this.addVariableWithElimination(2);
			Polynomial<T> t = ring.getVar(1);
			Polynomial<T> s = ring.getVar(2);
			Polynomial<T> sminust = ring.subtract(s, t);
			Polynomial<T> tsquare = ring.multiply(t, t);
			Polynomial<T> ssquare = ring.multiply(s, s);
			Polynomial<T> st = ring.multiply(s, t);
			generators.add(tsquare);
			generators.add(ssquare);
			generators.add(st);
			for (Polynomial<T> f : this.basis) {
				generators.add(ring.multiply(t, ring.getEmbeddingShift(f, 2)));
			}
			for (Polynomial<T> g : ideal.basis) {
				generators.add(ring.multiply(sminust, ring.getEmbeddingShift(g, 2)));
			}
			Ideal intersectionIdeal = ring.getIdeal(generators);
			List<Polynomial<T>> intersectionGenerators = new ArrayList<Polynomial<T>>();
			List<T> values = new ArrayList<T>();
			values.add(null);
			values.add(field.one());
			for (int i = 0; i < numvars; i++)
				values.add(null);
			for (Polynomial<T> b : intersectionIdeal.getBasis()) {
				if (b.getLeadingTerm().getExponents()[0] == 0 && b.getLeadingTerm().getExponents()[1] == 1)
					intersectionGenerators.add(PolynomialRing.this.getEmbeddingShift(b.evaluatePartially(values), -2));
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
		/*public int degree() {
			if (this.basis.isEmpty())
				return PolynomialRing.this.numvars + 1;
			int degree = 1;
			System.out.println("set: " + Arrays.toString(this.set));
			basisloop:
			for (int i = 0; i < this.basis.size(); i++) {
				int sum = 0;
				boolean nontrivial = false;
				for (int j = 0; j < this.leadingMonomials[i].length; j++) {
					System.out.print(this.leadingMonomials[i][j]);
					if (this.leadingMonomials[i][j] == 0)
						continue;
					if (this.set[j] != 0)
						continue basisloop;
					nontrivial = true;
					sum += this.leadingMonomials[i][j];
				}
				System.out.println();
				if (nontrivial)
					degree *= sum;
				System.out.println(sum);
			}
			return degree;
		}*/
		public String toString() {
			return this.basis.toString();
		}
	}
}
