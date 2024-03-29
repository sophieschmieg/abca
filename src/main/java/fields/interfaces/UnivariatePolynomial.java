package fields.interfaces;

import java.util.List;
import java.util.SortedSet;

import fields.local.Value;
import fields.polynomials.Monomial;

public interface UnivariatePolynomial<T extends Element<T>>  extends Polynomial<T> {

	UnivariatePolynomialRing<T> getPolynomialRing();

	Ring<T> getRing();

	int numberOfVariables();

	int degree();

	Value order();

	T univariateCoefficient(int degree);

	T coefficient(Monomial m);
	
	List<T> coefficients();

	T leadingCoefficient();

	Monomial leadingMonomial();

	SortedSet<Monomial> monomials();

	int compareTo(Polynomial<T> o);

	String toString();

	String toString(boolean ascending);

	String toString(String variable, boolean ascending);

	String toString(String[] variables, boolean ascending);

}