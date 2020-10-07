package fields.polynomials;

import java.util.List;
import java.util.SortedSet;

import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;

public class UnivariatePolynomial<T extends Element<T>> extends AbstractElement<Polynomial<T>> implements Polynomial<T> {
	private Ring<T> baseRing;
	private UnivariatePolynomialRing<T> polynomialRing;
	private List<T> coefficients;
	private List<T> traces;

	UnivariatePolynomial(UnivariatePolynomialRing<T> polynomialRing, List<T> coefficients) {
		this.polynomialRing = polynomialRing;
		this.baseRing = this.polynomialRing.getRing();
		int size = coefficients.size();
		while (size > 0 && coefficients.get(size - 1).equals(baseRing.zero())) {
			size--;
		}
		this.coefficients = coefficients.subList(0, size);
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
		return coefficients.size() - 1;
	}

	public T univariateCoefficient(int degree) {
		if (degree > degree()) {
			return this.baseRing.zero();
		}
		return this.coefficients.get(degree);
	}

	@Override
	public T coefficient(Monomial m) {
		return univariateCoefficient(m.exponents()[0]);
	}
	
	List<T> coefficients() {
		return this.coefficients;
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

	T trace(int degree) {
		if (degree < 1) {
			throw new ArithmeticException("traces are only defined for positive integers");
		}
		if (traces.size() >= degree) {
			return traces.get(degree - 1);
		}
		T trace = this.baseRing.multiply(-degree, univariateCoefficient(degree));
		for (int k = 1; k < degree; k++) {
			trace = this.baseRing.subtract(trace, this.baseRing.multiply(trace(k), univariateCoefficient(degree - k)));
		}
		if (this.traces.size() != degree - 1) {
			throw new ArithmeticException("Trace recursion went wrong");
		}
		this.traces.add(trace);
		return trace;
	}

	@Override
	public int compareTo(Polynomial<T> o) {
		if (!(o instanceof UnivariatePolynomial<?>)) {
			throw new ArithmeticException("Comparing against non univariate polynomial");
		}
		UnivariatePolynomial<T> p = (UnivariatePolynomial<T>) o;
		for (int i = degree(); i >= 0; i--) {
			int cmp = univariateCoefficient(i).compareTo((T) p.univariateCoefficient(i));
			if (cmp != 0) {
				return cmp;
			}
		}
		return 0;
	}

	@Override
	public String toString() {
		if (this.coefficients.isEmpty())
			return "0";
		StringBuffer buf = new StringBuffer();
		boolean first = true;
		for (int i = degree(); i >= 0; i--) {
			if (this.univariateCoefficient(i).equals(baseRing.zero())) {
				continue;
			}
			if (first) {
				first = false;
			} else {
				buf.append(" + ");
			}
			if (i == 0) {
				buf.append(this.univariateCoefficient(i).toString());
			} else {
				if (this.univariateCoefficient(i).equals(this.baseRing.one())) {
					buf.append("X");
				} else {
					buf.append(this.univariateCoefficient(i).toString() + "*X");
				}
				if (i > 1) {
					buf.append("^" + i);
				}
			}
		}
		return buf.toString();
	}

}
