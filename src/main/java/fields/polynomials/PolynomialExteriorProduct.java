package fields.polynomials;

import java.math.BigInteger;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.helper.AbstractModule;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.DifferentialForms.DifferentialForm;
import fields.vectors.Vector;

public class DifferentialForms<T extends Element<T>> extends AbstractModule<Polynomial<T>, DifferentialForm<T>> {
	private static class PrimitiveDifferentialForm<T extends Element<T>>
			implements Comparable<PrimitiveDifferentialForm<T>> {
		private Set<Integer> variables;
		private String[] variableNames;

		private PrimitiveDifferentialForm(Set<Integer> variables, String[] variableNames) {
			this.variables = variables;
			this.variableNames = variableNames;
		}

		@Override
		public String toString() {
			StringBuilder build = new StringBuilder();
			boolean first = true;
			for (int i : variables) {
				if (first) {
					first = false;
				} else {
					build.append(" ^ ");
				}
				build.append("d");
				build.append(variableNames[i - 1]);
			}
			return build.toString();
		}

		@Override
		public int compareTo(PrimitiveDifferentialForm<T> o) {
			for (int i = 0; i < variableNames.length; i++) {
				boolean hasI = variables.contains(i + 1);
				boolean otherHasI = o.variables.contains(i + 1);
				if (hasI && !otherHasI) {
					return 1;
				}
				if (!hasI && otherHasI) {
					return -1;
				}
			}
			return 0;
		}
	}

	public static class DifferentialForm<T extends Element<T>> extends AbstractElement<DifferentialForm<T>> {
		private List<Polynomial<T>> coefficients;
		private String[] variableNames;

		private DifferentialForm(List<Polynomial<T>> coefficients, String[] variableNames) {
			this.coefficients = coefficients;
			this.variableNames = variableNames;
		}

		@Override
		public String toString() {
			StringBuilder build = new StringBuilder();
			for (int i = 0; i < coefficients.size(); i++) {
				if (i != 0) {
					build.append(" + ");
				}
				build.append("(");
				build.append(coefficients.get(i));
				build.append(")*d");
				build.append(variableNames[i]);
			}
			return build.toString();
		}

		@Override
		public int compareTo(DifferentialForm<T> o) {
			for (int i = 0; i < coefficients.size(); i++) {
				int cmp = coefficients.get(i).compareTo(o.coefficients.get(i));
				if (cmp != 0) {
					return cmp;
				}
			}
			return 0;
		}
	}

	private PolynomialRing<T> polynomialRing;
	private int order;

	public DifferentialForms(PolynomialRing<T> polynomialRing, int order) {
		this.polynomialRing = polynomialRing;
	}

	@Override
	public PolynomialRing<T> getRing() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public DifferentialForm<T> zero() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public DifferentialForm<T> add(DifferentialForm<T> s1, DifferentialForm<T> s2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public DifferentialForm<T> negative(DifferentialForm<T> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public DifferentialForm<T> scalarMultiply(Polynomial<T> t, DifferentialForm<T> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isFree() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Ideal<Polynomial<T>> annihilator() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isLinearIndependent(List<DifferentialForm<T>> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isGeneratingModule(List<DifferentialForm<T>> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public List<List<Polynomial<T>>> nonTrivialCombinations(List<DifferentialForm<T>> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<DifferentialForm<T>> getModuleGenerators() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Vector<Polynomial<T>> asVector(DifferentialForm<T> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Exactness exactness() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public DifferentialForm<T> getRandomElement() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isFinite() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterator<DifferentialForm<T>> iterator() {
		// TODO Auto-generated method stub
		return null;
	}
}
