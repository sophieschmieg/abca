package fields.interfaces;

import java.util.SortedSet;

import fields.local.Value;
import fields.polynomials.Monomial;

public interface Polynomial<T extends Element<T>> extends Element<Polynomial<T>> {
	public int numberOfVariables();

	public PolynomialRing<T> getPolynomialRing();

	public Ring<T> getRing();

	public int degree();

	public int degree(int variable);

	public Value order();
	
	public Value order(int variable);

	public T coefficient(Monomial m);

	public T leadingCoefficient();

	public Monomial leadingMonomial();

	public SortedSet<Monomial> monomials();

	public String toString(String[] variables, boolean ascending);
}
