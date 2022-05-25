package fields.polynomials;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.SortedSet;

import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.local.Value;
import fields.vectors.Vector;

public class GenericUnivariatePolynomial<T extends Element<T>> extends AbstractElement<Polynomial<T>>
		implements Polynomial<T>, UnivariatePolynomial<T> {
	private Ring<T> baseRing;
	private GenericUnivariatePolynomialRing<T> polynomialRing;
	private List<T> coefficients;
	private List<T> traces;
	private T content;
	private Value order;

	GenericUnivariatePolynomial(GenericUnivariatePolynomialRing<T> polynomialRing, List<T> coefficients) {
		this.polynomialRing = polynomialRing;
		this.baseRing = this.polynomialRing.getRing();
		int degree = coefficients.size() - 1;
		while (degree >= 0 && coefficients.get(degree).equals(baseRing.zero())) {
			degree--;
		}
		this.coefficients = coefficients.subList(0, degree + 1);
		this.traces = new ArrayList<>();
		this.content = null;
		this.order = null;
	}

	@Override
	public PolynomialRing<T> getPolynomialRing() {
		return polynomialRing;
	}

	@Override
	public Ring<T> getRing() {
		return baseRing;
	}

	@Override
	public int numberOfVariables() {
		return 1;
	}

	@Override
	public int degree() {
		return coefficients().size() - 1;
	}

	@Override
	public int degree(int variable) {
		return degree();
	}

	@Override
	public Value order() {
		if (order == null) {
			order = Value.INFINITY;
			for (int i = 0; i <= degree(); i++) {
				if (!univariateCoefficient(i).equals(baseRing.zero())) {
					order = new Value(i);
					break;
				}
			}
		}
		return order;
	}

	@Override
	public Value order(int variable) {
		return order();
	}

	@Override
	public T univariateCoefficient(int degree) {
		if (degree > degree()) {
			return this.baseRing.zero();
		}
		return coefficients.get(degree);
	}

	@Override
	public T coefficient(Monomial m) {
		return univariateCoefficient(m.exponents()[0]);
	}

	List<T> coefficients() {
		return Collections.unmodifiableList(this.coefficients);
	}

	@Override
	public T leadingCoefficient() {
		return univariateCoefficient(Math.max(0, degree()));
	}

	@Override
	public Monomial leadingMonomial() {
		return polynomialRing.getMonomial(new int[] { Math.max(0, degree()) });
	}

	@Override
	public SortedSet<Monomial> monomials() {
		return polynomialRing.monomials(degree());
	}

	T content() {
		return content;
	}

	void setContent(T content) {
		this.content = content;
	}

	Vector<T> trace(int degree) {
		if (degree < 1) {
			throw new ArithmeticException("traces are only defined for positive integers");
		}
		if (traces.size() >= degree) {
			return new Vector<>(traces.subList(0, degree));
		}
		T trace = baseRing.multiply(-degree, univariateCoefficient(degree));
		for (int k = 1; k < degree; k++) {
			trace = baseRing.subtract(trace, baseRing.multiply(trace(k).get(k), univariateCoefficient(degree - k)));
		}
		if (this.traces.size() != degree - 1) {
			throw new ArithmeticException("Trace recursion went wrong");
		}
		this.traces.add(trace);
		return trace(degree);
	}

	@Override
	public int compareTo(Polynomial<T> o) {
		if (!(o instanceof GenericUnivariatePolynomial<?>)) {
			throw new ArithmeticException("Comparing against non univariate polynomial");
		}
		UnivariatePolynomial<T> p = (UnivariatePolynomial<T>) o;
		int degree = Math.max(degree(), o.degree());
		for (int i = degree; i >= 0; i--) {
			int cmp = univariateCoefficient(i).compareTo(p.univariateCoefficient(i));
			if (cmp != 0) {
				return cmp;
			}
		}
		return 0;
	}

	@Override
	public String toString() {
		return toString(polynomialRing.getVariableName(), false);
	}
	
	@Override
	public String toString(boolean ascending) {
		return toString(polynomialRing.getVariableName(), ascending);
	}

	@Override
	public String toString(String variable, boolean ascending) {
		if (this.coefficients.isEmpty())
			return "0";
		StringBuffer buf = new StringBuffer();
		boolean first = true;
		for (int index = 0; index <= degree(); index++) {
			int i = index;
			if (!ascending) {
				i = degree() - index;
			}
			if (this.univariateCoefficient(i).equals(baseRing.zero())) {
				continue;
			}
			if (first) {
				first = false;
			} else {
				buf.append(" + ");
			}
			String coefficient = this.univariateCoefficient(i).toString();
			if (coefficient.contains(" ")) {
				coefficient = "(" + coefficient + ")";
			}
			if (i == 0) {
				buf.append(coefficient);
			} else {
				if (this.univariateCoefficient(i).equals(this.baseRing.one())) {
					buf.append(variable);
				} else {
					buf.append(coefficient + "*" + variable);
				}
				if (i > 1) {
					buf.append("^" + i);
				}
			}
		}
		return buf.toString();
	}

	@Override
	public String toString(String[] variables, boolean ascending) {
		return toString(variables[0], ascending);
	}
}
