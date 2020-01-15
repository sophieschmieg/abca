package fields;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

public class Polynomial<T extends Element> implements Element {
	private PolynomialRing<T> ring;
	private Field<T> field;
	private int numvars;
	private int degree;
	private SortedMap<Monomial, T> coeff;
	private Comparator<Polynomial.Monomial> comparator;
	public final static Comparator<Polynomial.Monomial> LEX = new LexographicalOrder();
	public final static Comparator<Polynomial.Monomial> GRLEX = new GradedLexographicalOrder();
	public final static Comparator<Polynomial.Monomial> GREVLEX = new GradedReversedLexographicalOrder();
	
	public Polynomial(Map<Monomial,T> coeff, PolynomialRing<T> ring) {
		this.field = ring.getField();
		this.ring = ring;
		this.comparator = ring.getComparator();
		this.coeff = new TreeMap<Monomial, T>(this.comparator);
		this.degree = -1;
		this.numvars = this.ring.getNumVars();
		for (Monomial m : coeff.keySet()) {
			if (this.numvars != m.numvars)
				throw new ArithmeticException("wrong number of vars");
			if (!this.field.zero().equals(coeff.get(m))) {
				this.coeff.put(m, coeff.get(m));
				if (m.degree() > this.degree)
					this.degree = m.degree();
			}
		}
	}
	public int getDegree() {
		return this.degree;
	}
	public int getNumVars() {
		return this.numvars;
	}
	public boolean isHomogeneous() {
		for (Monomial m : this.coeff.keySet())
			if (m.degree() != this.degree)
				return false;
		return true;
	}
	public Map<Monomial, T> getCoefficients() {
		return Collections.unmodifiableMap(this.coeff);
	}
	public T getCoefficient(Monomial m) {
		if (this.coeff.containsKey(m))
			return this.coeff.get(m);
		return this.field.zero();
	}
	public T getLeadingCoefficient() {
		if (this.coeff.isEmpty())
			return this.field.zero();
		return this.getCoefficient(this.coeff.lastKey());
	}
	public Polynomial.Monomial getLeadingTerm() {
		if (this.coeff.isEmpty())
			return null;
		return this.coeff.lastKey();
	}
	public Polynomial<T> normalize() {
		T lc = this.getLeadingCoefficient();
		Map<Polynomial.Monomial,T> coeff = new TreeMap<Polynomial.Monomial, T>(this.comparator);
		for (Polynomial.Monomial m : this.coeff.keySet())
			coeff.put(m, this.field.divide(getCoefficient(m), lc));
		return new Polynomial<T>(coeff, this.ring);
	}
	@Override 
	public String toString() {
		if (this.coeff.isEmpty())
			return "0";
		StringBuffer buf = new StringBuffer();
		boolean first = true;
		SortedSet<Monomial> reversed = new TreeSet<Polynomial.Monomial>(Collections.reverseOrder(this.comparator));
		reversed.addAll(this.coeff.keySet());
		for (Monomial m : reversed) {
			if (first)
				first = false;
			else
				buf.append(" + ");
			if (m.degree() == 0)
				buf.append(this.coeff.get(m).toString());
			else if (this.coeff.get(m).equals(this.field.one()))
				buf.append(m.toString());
			else
				buf.append(this.coeff.get(m).toString() + "*" + m.toString());
		}
		return buf.toString();
	}
	@Override
	public boolean equals(Object O) {
		if (!(O instanceof Polynomial))
			return false;
		@SuppressWarnings("unchecked")
		Polynomial<T> p = (Polynomial<T>)O;
		for (Monomial m : this.coeff.keySet()) {
			if (!p.coeff.containsKey(m))
				return false;
			if (!this.coeff.get(m).equals(p.coeff.get(m)))
				return false;
		}
		return true;
	}
	
	public T evaluate(List<T> values) {
		Polynomial<T> result = this.evaluatePartially(values);
		if (result.getDegree() > 0)
			throw new ArithmeticException("not enough values");
		Monomial m = new Monomial(new int[values.size()]);
		if (result.getCoefficients().containsKey(m))
			return result.getCoefficients().get(m);
		else
			return this.field.zero();
	}

	public Polynomial<T> dehom(int coord) {
		if (coord < 1 || coord > this.numvars || !this.isHomogeneous())
			throw new ArithmeticException("Not possible");
		List<T> values = new ArrayList<T>();
		for (int i = 0; i < this.numvars; i++)
			if (i == coord - 1)
				values.add(this.field.one());
			else
				values.add(null);
		return this.evaluatePartially(values);
	}
	public Polynomial<T> evaluatePartially(List<T> values) {
		if (values.size() != this.numvars)
			throw new ArithmeticException("wrong number of vars");
		Map<Polynomial.Monomial, T> coeff = new TreeMap<Polynomial.Monomial, T>(this.comparator);
		for (Monomial m : this.getCoefficients().keySet()) {
			Monomial mp = this.dehomMonomialPartially(m, values);
			if (!coeff.containsKey(mp))
				coeff.put(mp, this.field.zero());
			coeff.put(mp, this.field.add(coeff.get(mp), this.field.multiply(this.getCoefficients().get(m), this.evaluateMonomialPartially(m, values))));
		}
		return new Polynomial<T>(coeff, this.ring);
	}
	private T evaluateMonomialPartially(Monomial m, List<T> values) {
		T result = this.field.one();
		for (int i = 0; i < this.numvars; i++)
			if (values.get(i) != null)
				result = this.field.multiply(result, this.field.power(values.get(i), m.getExponents()[i]));
		return result;
	}
	private Polynomial.Monomial dehomMonomialPartially(Monomial m, List<T> values) {
		int[] exponents = new int[numvars];
		for (int i = 0; i < this.numvars; i++)
			if (values.get(i) == null)
				exponents[i] = m.getExponents()[i];
		return new Polynomial.Monomial(exponents);
	}
	public Polynomial<T> derivative(int coord) {
		if (coord < 1 || coord > this.numvars)
			throw new ArithmeticException("wrong number of vars");
		coord--;
		Map<Polynomial.Monomial, T> coeff = new TreeMap<Polynomial.Monomial, T>(this.comparator);
		coeff.put(new Monomial(new int[this.numvars]), this.field.zero());
		for (Monomial m : this.getCoefficients().keySet()) {
			int k = m.getExponents()[coord];
			if (k == 0)
				continue;
			int[] exponents = new int[this.numvars];
			for (int i = 0; i < this.numvars; i++)
				exponents[i] = m.getExponents()[i];
			exponents[coord]--;
			Monomial mp = new Monomial(exponents);
			if (!coeff.containsKey(mp))
				coeff.put(mp, this.field.zero());
			coeff.put(mp, this.field.add(coeff.get(mp), this.field.multiply(k, this.getCoefficients().get(m))));
		}
		return new Polynomial<T>(coeff, this.ring);
	}

	@Override
	public int compareTo(Element e) {
		@SuppressWarnings("unchecked")
		Polynomial<T> p = (Polynomial<T>)e;
		Set<Monomial> set = new TreeSet<Monomial>(Collections.reverseOrder(this.comparator));
		set.addAll(this.coeff.keySet());
		set.addAll(p.coeff.keySet());
		for (Monomial m : set) {
			int cmp = this.getCoefficient(m).compareTo(p.getCoefficient(m));
			if (cmp != 0)
				return cmp;
		}
		return 0;
	}
	public PolynomialRing<T> getRing() {
		return this.ring;
	}
	public static class Monomial {
		private int numvars;
		private int[] exponents;
		public Monomial(int[] exponents) {
			this.numvars = exponents.length;
			this.exponents = Arrays.copyOf(exponents, this.numvars);
			for (int e : this.exponents)
				if (e < 0)
					throw new ArithmeticException("negative exponent");
		}
		@Override
		public boolean equals(Object O) {
			if (!(O instanceof Monomial))
				return false;
			Monomial m = (Monomial)O;
			if (this.numvars != m.numvars)
				return false;
			for (int i = 0; i < this.numvars; i++)
				if (this.exponents[i] != m.exponents[i])
					return false;
			return true;
		}
		
		public int degree() {
			int d = 0;
			for (int e : this.exponents)
				d += e;
			return d;
		}
		public static Monomial multiply(Monomial t1, Monomial t2) {
			int numvars = t1.numvars;
			if (numvars != t2.numvars)
				throw new ArithmeticException("wrong number of variables");
			int[] exponents = new int[numvars];
			for (int i = 0; i < numvars; i++)
				exponents[i] = t1.exponents[i] + t2.exponents[i];
			return new Monomial(exponents);
		}
		public int[] getExponents() {
			return this.exponents;
		}
		@Override
		public String toString() {
			StringBuffer buf = new StringBuffer();
			String[] vars;
			if (this.numvars == 1)
				vars = new String[] {"X"};
			else if (this.numvars == 2)
				vars = new String[] {"X", "Y"};
			else if (this.numvars == 3)
				vars = new String[] {"X", "Y", "Z"};
			else {
				vars = new String[this.numvars];
				for (int i = 0; i < this.numvars; i++)
					vars[i] = "X_" + i;
			}
			for (int i = 0; i < this.numvars; i++) {
				if (this.exponents[i] > 1)
					buf.append(vars[i] + "^" + this.exponents[i]);
				else if (this.exponents[i] == 1)
					buf.append(vars[i]);
			}
			return buf.toString();
		}
	}
	public static class LexographicalOrder implements Comparator<Polynomial.Monomial> {

		@Override
		public int compare(Monomial a, Monomial b) {
			if (a.numvars != b.numvars)
				return a.numvars - b.numvars;
			for (int i = 0; i < a.numvars; i++) {
				if (a.exponents[i] != b.exponents[i])
					return a.exponents[i] - b.exponents[i];
			}
			return 0;
		}
		
	}
	public static class GradedLexographicalOrder implements Comparator<Polynomial.Monomial> {

		@Override
		public int compare(Monomial a, Monomial b) {
			if (a.numvars != b.numvars)
				return a.numvars - b.numvars;
			if (a.degree() != b.degree())
				return a.degree() - b.degree();
			for (int i = 0; i < a.numvars; i++) {
				if (a.exponents[i] != b.exponents[i])
					return a.exponents[i] - b.exponents[i];
			}
			return 0;
		}
	}
	public static class GradedReversedLexographicalOrder implements Comparator<Polynomial.Monomial> {

		@Override
		public int compare(Monomial a, Monomial b) {
			if (a.numvars != b.numvars)
				return a.numvars - b.numvars;
			if (a.degree() != b.degree())
				return a.degree() - b.degree();
			for (int i = 0; i < a.numvars; i++) {
				if (a.exponents[a.numvars - i - 1] != b.exponents[a.numvars - i - 1])
					return b.exponents[a.numvars - i - 1] - a.exponents[a.numvars - i - 1];
			}
			return 0;
		}
	}
	public static class EliminationOrder implements Comparator<Polynomial.Monomial> {
		private Comparator<Polynomial.Monomial> firstcomp;
		private Comparator<Polynomial.Monomial> secondcomp;
		private int blocklength;
		public EliminationOrder(Comparator<Monomial> firstcomp,
				Comparator<Monomial> secondcomp, int blocklength) {
			this.firstcomp = firstcomp;
			this.secondcomp = secondcomp;
			this.blocklength = blocklength;
		}
		@Override
		public int compare(Monomial a, Monomial b) {
			if (a.numvars != b.numvars)
				return a.numvars - b.numvars;
			if (a.numvars <= this.blocklength)
				return this.firstcomp.compare(a, b);
			int[] afirstblock = Arrays.copyOf(a.exponents, this.blocklength);
			int[] bfirstblock = Arrays.copyOf(b.exponents, this.blocklength);
			int cmp = this.firstcomp.compare(new Monomial(afirstblock), new Monomial(bfirstblock));
			if (cmp != 0)
				return cmp;
			int[] asecondblock = Arrays.copyOfRange(a.exponents, this.blocklength, a.numvars);
			int[] bsecondblock = Arrays.copyOfRange(b.exponents, this.blocklength, b.numvars);
			return this.secondcomp.compare(new Monomial(asecondblock), new Monomial(bsecondblock));
		}
	}
}
