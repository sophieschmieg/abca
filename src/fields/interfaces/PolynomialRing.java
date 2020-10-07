package fields.interfaces;

import java.util.Comparator;
import java.util.List;
import java.util.Map;

import fields.polynomials.Monomial;
import varieties.AffinePoint;

public interface PolynomialRing<T extends Element<T>> extends Algebra<T, Polynomial<T>> {
	public int numberOfVariables();
	public Comparator<Monomial> getComparator();
	public Polynomial<T> getPolynomial(Map<Monomial, T> t);
	public Monomial getMonomial(int[] exponents);
	public Polynomial<T> getVar(int i);
	public Polynomial<T> getVarPower(int i, int power);
	public Polynomial<T> getEmbedding(T t, int[] exponents);
	public Polynomial<T> getLinear(List<T> coeff);
	public Polynomial<T> getLinear(@SuppressWarnings("unchecked") T... coeff);

	public PolynomialRing<T> addVariableWithElimination(int shift);
	public Polynomial<T> getEmbeddingWithElimination(Polynomial<T> t, int shift);
	public Polynomial<T> getEmbedding(Polynomial<T> t, int[] map);

	public Polynomial<T> multiply(T t1, Polynomial<T> t2);
	public Polynomial<T> multiply(T t1, Monomial m, Polynomial<T> t2);
	public Polynomial<T> divide(Polynomial<T> dividend, T divisor);
	public List<Polynomial<T>> generalQuotientAndRemainder(Polynomial<T> dividend, List<Polynomial<T>> divisors);
	public Polynomial<T> reduce(Polynomial<T> polynomial, List<Polynomial<T>> basis);
	public List<Polynomial<T>> buchberger(List<Polynomial<T>> generators);
	public List<Polynomial<T>> reduceBasis(List<Polynomial<T>> generators);
	public List<AffinePoint<T>> solve(List<Polynomial<T>> polynomials);
		
	public T evaluate(Polynomial<T> t, @SuppressWarnings("unchecked") T... ts);
	public T evaluate(Polynomial<T> t, List<T> ts);
	public Polynomial<T> partiallyEvaluate(Polynomial<T> t, List<T> ts);
	public Polynomial<T> substitute(Polynomial<T> t, List<Polynomial<T>> values);

	public boolean canBeNormalized(Polynomial<T> t);
	public Polynomial<T> normalize(Polynomial<T> t);
	public T content(Polynomial<T> t);
	public Polynomial<T> contentFree(Polynomial<T> t);

	public Polynomial<T> derivative(Polynomial<T> t, int variable);
	
	public boolean isHomogeneous(Polynomial<T> t);
	public Polynomial<T> homogenize(Polynomial<T> t);
	public Polynomial<T> dehomogenize(Polynomial<T> t, int coord);
		
	public Polynomial<T> round(Polynomial<T> t, int degree);
}
