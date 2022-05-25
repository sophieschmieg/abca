package fields.polynomials;

import java.util.Arrays;
import java.util.Comparator;

public class Monomial implements Comparable<Monomial> {
	private Comparator<Monomial> comparator;
	private int numvars;
	private int[] exponents;

	public final static Comparator<Monomial> LEX = new LexographicalOrder();
	public final static Comparator<Monomial> REVLEX = new ReverseLexographicalOrder();
	public final static Comparator<Monomial> GRLEX = new GradedLexographicalOrder();
	public final static Comparator<Monomial> GREVLEX = new GradedReversedLexographicalOrder();

	public Monomial(Comparator<Monomial> comparator, int[] exponents) {
		this.comparator = comparator;
		this.numvars = exponents.length;
		this.exponents = Arrays.copyOf(exponents, this.numvars);
		for (int e : this.exponents)
			if (e < 0)
				throw new ArithmeticException("negative exponent");
	}

	@Override
	public int compareTo(Monomial o) {
		return this.comparator.compare(this, o);
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Monomial))
			return false;
		return this.compareTo((Monomial) o) == 0;
	}

	public int degree() {
		int d = 0;
		for (int e : this.exponents) {
			d += e;
		}
		return d;
	}

	Monomial homogenizeMonomial(int degree) {
		int[] exponents = Arrays.copyOf(this.exponents, this.numvars + 1);
		exponents[this.numvars] = degree - this.degree();
		return new Monomial(comparator, exponents);
	}

	Monomial homogenizeMonomial(int degree, int coord) {
		int[] exponents = Arrays.copyOf(this.exponents, this.numvars);
		if (exponents[coord - 1] != 0) {
			throw new ArithmeticException("not dehomogenized in coordinate " + coord);
		}
		exponents[coord - 1] = degree - this.degree();
		return new Monomial(comparator, exponents);
	}

	public Monomial multiply(Monomial m) {
		if (numvars != m.numvars) {
			throw new ArithmeticException("wrong number of variables");
		}
		int[] exponents = new int[numvars];
		for (int i = 0; i < numvars; i++) {
			exponents[i] = this.exponents[i] + m.exponents[i];
		}
		return new Monomial(comparator, exponents);
	}

	public Monomial divide(Monomial m) {
		int[] dividedexponents = new int[this.numvars];
		for (int i = 0; i < this.numvars; i++) {
			if (exponents()[i] < m.exponents()[i]) {
				return null;
			}
			dividedexponents[i] = exponents()[i] - m.exponents()[i];
		}
		return new Monomial(comparator, dividedexponents);
	}

	public Monomial lcm(Monomial m) {
		int[] exponents = new int[this.numvars];
		for (int i = 0; i < this.numvars; i++) {
			exponents[i] = Math.max(exponents()[i], m.exponents()[i]);
		}
		return new Monomial(comparator, exponents);
	}

	public int[] exponents() {
		return this.exponents;
	}

	public String[] defaultVariableNames() {
		String[] vars;
		if (this.numvars == 1)
			vars = new String[] { "X" };
		else if (this.numvars == 2)
			vars = new String[] { "X", "Y" };
		else if (this.numvars == 3)
			vars = new String[] { "X", "Y", "Z" };
		else {
			vars = new String[this.numvars];
			for (int i = 0; i < this.numvars; i++)
				vars[i] = "X_" + i;
		}
		return vars;
	}

	@Override
	public String toString() {
		return toString(defaultVariableNames());
	}

	public String toString(String[] variables) {
		StringBuffer buf = new StringBuffer();
		boolean first = true;
		for (int i = 0; i < this.numvars; i++) {
			if (first && exponents[i] != 0) {
				first=false;
			} else if (exponents[i] != 0) {
				buf.append("*");
			}
			if (this.exponents[i] > 1)
				buf.append(variables[i] + "^" + this.exponents[i]);
			else if (this.exponents[i] == 1)
				buf.append(variables[i]);
		}
		return buf.toString();
	}
	
	public static class LexographicalOrder implements Comparator<Monomial> {

		@Override
		public int compare(Monomial a, Monomial b) {
			if (a.numvars != b.numvars)
				throw new ArithmeticException("number of variables not aligned");
			for (int i = 0; i < a.numvars; i++) {
				if (a.exponents[i] != b.exponents[i])
					return a.exponents[i] - b.exponents[i];
			}
			return 0;
		}

	}

	public static class ReverseLexographicalOrder implements Comparator<Monomial> {

		@Override
		public int compare(Monomial a, Monomial b) {
			if (a.numvars != b.numvars)
				throw new ArithmeticException("number of variables not aligned");
			for (int i = 0; i < a.numvars; i++) {
				if (a.exponents[a.numvars - i - 1] != b.exponents[a.numvars - i - 1])
					return a.exponents[a.numvars - i - 1] - b.exponents[a.numvars - i - 1];
			}
			return 0;
		}

	}

	public static class GradedLexographicalOrder implements Comparator<Monomial> {

		@Override
		public int compare(Monomial a, Monomial b) {
			if (a.numvars != b.numvars)
				throw new ArithmeticException("number of variables not aligned");
			if (a.degree() != b.degree())
				return a.degree() - b.degree();
			for (int i = 0; i < a.numvars; i++) {
				if (a.exponents[i] != b.exponents[i])
					return a.exponents[i] - b.exponents[i];
			}
			return 0;
		}
	}

	public static class GradedReversedLexographicalOrder implements Comparator<Monomial> {

		@Override
		public int compare(Monomial a, Monomial b) {
			if (a.numvars != b.numvars)
				throw new ArithmeticException("number of variables not aligned");
			if (a.degree() != b.degree())
				return a.degree() - b.degree();
			for (int i = 0; i < a.numvars; i++) {
				if (a.exponents[a.numvars - i - 1] != b.exponents[a.numvars - i - 1])
					return b.exponents[a.numvars - i - 1] - a.exponents[a.numvars - i - 1];
			}
			return 0;
		}
	}

	public static class EliminateVariableOrder implements Comparator<Monomial> {
		private int variable;

		public EliminateVariableOrder(int variable) {
			this.variable = variable - 1;
		}

		@Override
		public int compare(Monomial a, Monomial b) {
			if (a.numvars != b.numvars) {
				throw new ArithmeticException("number of variables not aligned");
			}
			if (a.exponents[variable] != b.exponents[variable]) {
				return a.exponents[variable] - b.exponents[variable];
			}
			if (a.degree() != b.degree()) {
				return a.degree() - b.degree();
			}
			for (int i = a.numvars - 1; i >= 0; i--) {
				if (a.exponents[i] != b.exponents[i]) {
					return a.exponents[i] - b.exponents[i];
				}
			}
			return 0;
		}

	}

	public static class InvertedEliminateVariableOrder implements Comparator<Monomial> {
		private int variable;

		public InvertedEliminateVariableOrder(int variable) {
			this.variable = variable - 1;
		}

		@Override
		public int compare(Monomial a, Monomial b) {
			if (a.numvars != b.numvars) {
				throw new ArithmeticException("number of variables not aligned");
			}
			if (a.degree() - a.exponents[variable] != b.degree() - b.exponents[variable]) {
				return a.degree() - a.exponents[variable] - b.degree() + b.exponents[variable];
			}
			for (int i = a.numvars - 1; i >= 0; i--) {
				if (i == variable) {
					continue;
				}
				if (a.exponents[i] != b.exponents[i]) {
					return a.exponents[i] - b.exponents[i];
				}
			}
			if (a.exponents[variable] != b.exponents[variable]) {
				return a.exponents[variable] - b.exponents[variable];
			}
			return 0;
		}

	}

	public static class EliminationOrder implements Comparator<Monomial> {
		private Comparator<Monomial> firstcomp;
		private Comparator<Monomial> secondcomp;
		private int blocklength;

		public EliminationOrder(Comparator<Monomial> firstcomp, Comparator<Monomial> secondcomp, int blocklength) {
			this.firstcomp = firstcomp;
			this.secondcomp = secondcomp;
			this.blocklength = blocklength;
		}

		@Override
		public int compare(Monomial a, Monomial b) {
			if (a.numvars != b.numvars)
				throw new ArithmeticException("number of variables not aligned");
			if (a.numvars <= this.blocklength)
				return this.firstcomp.compare(a, b);
			int[] afirstblock = Arrays.copyOf(a.exponents, this.blocklength);
			int[] bfirstblock = Arrays.copyOf(b.exponents, this.blocklength);
			int cmp = this.firstcomp.compare(new Monomial(this.firstcomp, afirstblock),
					new Monomial(this.firstcomp, bfirstblock));
			if (cmp != 0)
				return cmp;
			int[] asecondblock = Arrays.copyOfRange(a.exponents, this.blocklength, a.numvars);
			int[] bsecondblock = Arrays.copyOfRange(b.exponents, this.blocklength, b.numvars);
			return this.secondcomp.compare(new Monomial(this.secondcomp, asecondblock),
					new Monomial(this.secondcomp, bsecondblock));
		}
	}

}
