package fields.polynomials;

import java.util.Collections;
import java.util.Comparator;
import java.util.Map;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring;
import fields.local.Value;
import util.SortedKeySet;

public class MultivariatePolynomial<T extends Element<T>> extends AbstractElement<Polynomial<T>>
		implements Polynomial<T> {
	private MultivariatePolynomialRing<T> polynomialRing;
	private Ring<T> baseRing;
	private int numvars;
	private int degree;
	private SortedMap<Monomial, T> coeff;

	MultivariatePolynomial(Map<Monomial, T> coeff, MultivariatePolynomialRing<T> ring) {
		this.baseRing = ring.getRing();
		this.polynomialRing = ring;
		this.coeff = new TreeMap<>();
		this.degree = -1;
		this.numvars = this.polynomialRing.numberOfVariables();
		for (Monomial m : coeff.keySet()) {
			if (this.numvars != m.exponents().length) {
				throw new ArithmeticException("wrong number of vars");
			}
			if (!this.baseRing.zero().equals(coeff.get(m))) {
				this.coeff.put(m, coeff.get(m));
				if (m.degree() > this.degree) {
					this.degree = m.degree();
				}
			}
		}
	}

	@Override
	public int degree() {
		return this.degree;
	}

	@Override
	public int degree(int variable) {
		int degree = -1;
		for (Monomial m : coeff.keySet()) {
			if (degree < m.exponents()[variable - 1]) {
				degree = m.exponents()[variable - 1];
			}
		}
		return degree;
	}

	@Override
	public Value order() {
		Value min = Value.INFINITY;
		for (Monomial m : coeff.keySet()) {
			Value v = new Value(m.degree());
			if (min.compareTo(v) > 0) {
				min = v;
			}
		}
		return min;
	}

	@Override
	public Value order(int variable) {
		Value min = Value.INFINITY;
		for (Monomial m : coeff.keySet()) {
			Value v = new Value(m.exponents()[variable - 1]);
			if (min.compareTo(v) > 0) {
				min = v;
			}
		}
		return min;
	}

	@Override
	public int numberOfVariables() {
		return this.numvars;
	}

	@Override
	public SortedSet<Monomial> monomials() {
		return new SortedKeySet<Monomial, T>(this.coeff);
	}

	@Override
	public T coefficient(Monomial m) {
		if (this.coeff.containsKey(m))
			return this.coeff.get(m);
		return this.baseRing.zero();
	}

	@Override
	public T leadingCoefficient() {
		if (this.coeff.isEmpty())
			return this.baseRing.zero();
		return this.coefficient(this.coeff.lastKey());
	}

	public Monomial leadingMonomial() {
		if (this.coeff.isEmpty()) {
			return polynomialRing.getMonomial(new int[numvars]);
		}
		return this.coeff.lastKey();
	}

	@Override
	public String toString() {
		if (this.coeff.isEmpty()) {
			return "0";
		}
		return toString(polynomialRing.getVariableNames(), false);
	}

	public String toString(String[] variables, boolean ascending) {
		if (this.coeff.isEmpty()) {
			return "0";
		}
		StringBuffer buf = new StringBuffer();
		boolean first = true;
		SortedSet<Monomial> reversed = new TreeSet<Monomial>(
				ascending ? Comparator.naturalOrder() : Collections.reverseOrder());
		reversed.addAll(this.coeff.keySet());
		for (Monomial m : reversed) {
			if (first)
				first = false;
			else
				buf.append(" + ");
			if (m.degree() == 0)
				buf.append(this.coeff.get(m).toString());
			else if (this.coeff.get(m).equals(this.baseRing.one()))
				buf.append(m.toString(variables));
			else
				buf.append(this.coeff.get(m).toString() + "*" + m.toString(variables));
		}
		return buf.toString();
	}

	@Override
	public int compareTo(Polynomial<T> p) {
		SortedSet<Monomial> set = new TreeSet<Monomial>(Collections.reverseOrder());
		set.addAll(this.monomials());
		set.addAll(p.monomials());
		for (Monomial m : set) {
			int cmp = this.coefficient(m).compareTo(p.coefficient(m));
			if (cmp != 0)
				return cmp;
		}
		return 0;
	}

	public Ring<T> getRing() {
		return baseRing;
	}

	public MultivariatePolynomialRing<T> getPolynomialRing() {
		return this.polynomialRing;
	}
}
