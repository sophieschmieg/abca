package fields.interfaces;

import java.util.Comparator;
import java.util.List;
import java.util.Map;

import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.Vector;
import varieties.affine.AffinePoint;

public interface PolynomialRing<T extends Element<T>> extends Algebra<T, Polynomial<T>> {
	public int numberOfVariables();
	public Comparator<Monomial> getComparator();
	public Polynomial<T> getPolynomial(Map<Monomial, T> t);
	public Monomial getMonomial(int[] exponents);
	public Polynomial<T> getVar(int i);
	public Polynomial<T> getVarPower(int i, int power);
	public Polynomial<T> getEmbedding(T t, int[] exponents);
	public <S extends Element<S>> Polynomial<T> getEmbedding(Polynomial<S> t, MathMap<S, T> mapping);
	public Polynomial<T> getLinear(List<T> coeff);
	public Polynomial<T> getLinear(@SuppressWarnings("unchecked") T... coeff);
	public Polynomial<T> getRandomElement(int degree);

	public PolynomialRing<T> eliminateVariable();
	public UnivariatePolynomial<Polynomial<T>> asUnivariatePolynomial(Polynomial<T> t, int variable);
	public Polynomial<T> fromUnivariatePolynomial(UnivariatePolynomial<Polynomial<T>> t, int variable);
	public PolynomialRing<T> addVariableWithElimination(int shift);
	public Polynomial<T> getEmbeddingWithElimination(Polynomial<T> t, int shift);
	public Polynomial<T> getEmbedding(Polynomial<T> t, int[] map);
	public Polynomial<T> getEmbedding(Polynomial<T> t);

	public Polynomial<T> multiply(T t1, Polynomial<T> t2);
	public Polynomial<T> multiply(T t1, Monomial m, Polynomial<T> t2);
	public Polynomial<T> divideScalar(Polynomial<T> dividend, T divisor);
	public static class GeneralQuotientAndRemainderResult<T extends Element<T>> {
		private Polynomial<T> remainder;
		private List<Polynomial<T>> quotients;
		
		public GeneralQuotientAndRemainderResult(Polynomial<T> remainder, List<Polynomial<T>> quotients) {
			this.remainder = remainder;
			this.quotients = quotients;
		}
		public Polynomial<T> getRemainder() {
			return remainder;
		}
		public List<Polynomial<T>> getQuotients() {
			return quotients;
		}
	}
	public GeneralQuotientAndRemainderResult<T> generalQuotientAndRemainder(Polynomial<T> dividend, List<Polynomial<T>> divisors);
	public Polynomial<T> reduce(Polynomial<T> polynomial, List<Polynomial<T>> basis);
	public static class GroebnerBasis<T extends Element<T>> {
		private List<Polynomial<T>> basis;
		private List<List<Polynomial<T>>> expression;

		public GroebnerBasis(List<Polynomial<T>> basis, List<List<Polynomial<T>>> expression) {
			this.basis = basis;
			this.expression = expression;
		}

		public List<Polynomial<T>> getBasis() {
			return basis;
		}

		public List<List<Polynomial<T>>> getExpression() {
			return expression;
		}
	}
	public GroebnerBasis<T> buchberger(List<Polynomial<T>> generators);
	public GroebnerBasis<T> reduceBasis(GroebnerBasis<T> basis);
	
	public List<AffinePoint<T>> solve(List<Polynomial<T>> polynomials);
	
	public IdealResult<Polynomial<T>, PolynomialIdeal<T>> getIdealWithTransforms(List<Polynomial<T>> generators);
	public IdealResult<Polynomial<T>, PolynomialIdeal<T>> getIdealWithTransforms(@SuppressWarnings("unchecked") Polynomial<T>... generators);
	public PolynomialIdeal<T> getIdeal(List<Polynomial<T>> t);
	public PolynomialIdeal<T> getIdeal(@SuppressWarnings("unchecked") Polynomial<T>... t);
	public PolynomialIdeal<T> getUnitIdeal();
	public PolynomialIdeal<T> getZeroIdeal();
	public PolynomialIdeal<T> add(Ideal<Polynomial<T>> t1, Ideal<Polynomial<T>> t2);
	public PolynomialIdeal<T> intersect(Ideal<Polynomial<T>> t1, Ideal<Polynomial<T>> t2);
	
	public T evaluate(Polynomial<T> t, @SuppressWarnings("unchecked") T... ts);
	public T evaluate(Polynomial<T> t, List<T> ts);
	public T evaluate(Polynomial<T> t, Vector<T> ts);
	public T evaluate(Polynomial<T> t, AffinePoint<T> ts);
	public Polynomial<T> partiallyEvaluate(Polynomial<T> t, List<T> ts);
	public Polynomial<T> substitute(Polynomial<T> t, List<Polynomial<T>> values);
	
	public boolean canBeNormalized(Polynomial<T> t);
	public Polynomial<T> normalize(Polynomial<T> t);
	public T content(Polynomial<T> t);
	public Polynomial<T> contentFree(Polynomial<T> t);
	public FactorizationResult<Polynomial<T>, T> squareFreeFactorization(Polynomial<T> t);

	public Polynomial<T> derivative(Polynomial<T> t, int variable);
	
	public Polynomial<T> resultant(Polynomial<T> t1, Polynomial<T> t2, int variable);
	
	public boolean isHomogeneous(Polynomial<T> t);
	public Polynomial<T> homogenize(Polynomial<T> t);
	public Polynomial<T> homogenize(Polynomial<T> t, int coord);
	public Polynomial<T> dehomogenize(Polynomial<T> t, int coord);
			
	public Polynomial<T> round(Polynomial<T> t, int degree);
}
