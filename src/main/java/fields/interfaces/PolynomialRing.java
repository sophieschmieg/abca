package fields.interfaces;

import java.util.Comparator;
import java.util.List;
import java.util.Map;

import fields.polynomials.DifferentialForms;
import fields.polynomials.DifferentialForms.DifferentialForm;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.Matrix;
import fields.vectors.Vector;
import varieties.affine.AffinePoint;

public interface PolynomialRing<T extends Element<T>> extends Algebra<T, Polynomial<T>> {
	public int numberOfVariables();
	public PolynomialRing<T> withVariableNames(String[] variableNames);
	public void setVariableNames(String[] variableNames);
	public String[] getVariableNames();
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
	public ModuloMaximalIdealResult<Polynomial<T>, ?, PolynomialRing<T>, PolynomialIdeal<T>, ?> moduloMaximalIdeal(Ideal<Polynomial<T>> ideal);

	public PolynomialRing<T> eliminateVariable();
	public PolynomialRing<T> eliminateVariable(int variable);
	public UnivariatePolynomial<Polynomial<T>> asUnivariatePolynomial(Polynomial<T> t, int variable);
	public Polynomial<T> fromUnivariatePolynomial(UnivariatePolynomial<Polynomial<T>> t, int variable);
	public Polynomial<T> flattenPolynomial(Polynomial<Polynomial<T>> t);
	public Polynomial<Polynomial<T>> unflattenPolynomial(Polynomial<T> t, PolynomialRing<Polynomial<T>> polynomialRing, PolynomialRing<T> baseRing);
    public PolynomialRing<T> addVariableWithElimination(int shift);
	public Polynomial<T> getEmbeddingWithElimination(Polynomial<T> t, int shift);
	public Polynomial<T> getEmbedding(Polynomial<T> t, int[] map);
	public Polynomial<T> getEmbedding(Polynomial<T> t);

	public Polynomial<T> multiply(T t1, Polynomial<T> t2);
	public Polynomial<T> multiply(T t1, Polynomial<T> t2, Polynomial<T> t3);
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
		private List<Vector<Polynomial<T>>> syzygies;

		public GroebnerBasis(List<Polynomial<T>> basis, List<List<Polynomial<T>>> expression, List<Vector<Polynomial<T>>> syzygies) {
			this.basis = basis;
			this.expression = expression;
			this.syzygies = syzygies;
		}

		public List<Polynomial<T>> getBasis() {
			return basis;
		}

		public List<List<Polynomial<T>>> getExpression() {
			return expression;
		}
		
		public List<Vector<Polynomial<T>>> getSyzygies() {
			return syzygies;
		}
	}
	public GroebnerBasis<T> buchberger(List<Polynomial<T>> generators, boolean computeExpressionsAndSyzygies);
	public GroebnerBasis<T> reduceBasis(GroebnerBasis<T> basis, boolean computeExpressionsAndSyzygies);
	
	public List<AffinePoint<T>> solve(List<Polynomial<T>> polynomials);
	public List<AffinePoint<T>> solve(PolynomialIdeal<T> ideal);
	
	public IdealResult<Polynomial<T>, PolynomialIdeal<T>> getIdealWithTransforms(List<Polynomial<T>> generators);
	public IdealResult<Polynomial<T>, PolynomialIdeal<T>> getIdealWithTransforms(@SuppressWarnings("unchecked") Polynomial<T>... generators);
	public PolynomialIdeal<T> getIdeal(List<Polynomial<T>> t);
	public PolynomialIdeal<T> getIdeal(@SuppressWarnings("unchecked") Polynomial<T>... t);
	public PolynomialIdeal<T> getUnitIdeal();
	public PolynomialIdeal<T> getZeroIdeal();
	public PolynomialIdeal<T> add(Ideal<Polynomial<T>> t1, Ideal<Polynomial<T>> t2);
	public PolynomialIdeal<T> intersect(Ideal<Polynomial<T>> t1, Ideal<Polynomial<T>> t2);
	public PolynomialIdeal<T> radical(Ideal<Polynomial<T>> t);
	public PolynomialIdeal<T> homogenizeIdeal(Ideal<Polynomial<T>> ideal);
	public PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> primaryDecomposition(Ideal<Polynomial<T>> t);
	public PolynomialIdeal<T> getEmbeddingOfBaseIdeal(Ideal<T> t);
	public PolynomialIdeal<T> getEmbedding(Ideal<Polynomial<T>> t);
	public PolynomialIdeal<T> getEmbedding(Ideal<Polynomial<T>> t, int[] map);
	public <S extends Element<S>> PolynomialIdeal<T> getEmbedding(Ideal<Polynomial<S>> t, MathMap<S, T> map);
	public List <PolynomialIdeal<T>> maximalPrimeIdealChain();
	public List <PolynomialIdeal<T>> maximalPrimeIdealChain(Ideal<Polynomial<T>> start);
	public List <PolynomialIdeal<T>> maximalPrimeIdealChain(Ideal<Polynomial<T>> start,Ideal<Polynomial<T>> end);
	
	public T evaluate(Polynomial<T> t, @SuppressWarnings("unchecked") T... ts);
	public T evaluate(Polynomial<T> t, List<T> ts);
	public T evaluate(Polynomial<T> t, Vector<T> ts);
	public T evaluate(Polynomial<T> t, AffinePoint<T> ts);
	public Vector<T> evaluate(Vector<Polynomial<T>> t, Vector<T> ts);
	public Matrix<T> evaluate(Matrix<Polynomial<T>> t, Vector<T> ts);
	public Polynomial<T> partiallyEvaluate(Polynomial<T> t, List<T> ts);
	public Polynomial<T> substitute(Polynomial<T> t, List<Polynomial<T>> values);
	
	public boolean canBeNormalized(Polynomial<T> t);
	public Polynomial<T> normalize(Polynomial<T> t);
	public T content(Polynomial<T> t);
	public Polynomial<T> contentFree(Polynomial<T> t);
	public FactorizationResult<Polynomial<T>, T> squareFreeFactorization(Polynomial<T> t);
	public Polynomial<T> radical(Polynomial<T> t);

	public Polynomial<T> derivative(Polynomial<T> t, int variable);
	public Vector<Polynomial<T>> gradient(Polynomial<T> t);
	public Matrix<Polynomial<T>> jacobianMatrix(Vector<Polynomial<T>> t);
	public DifferentialForm<T> totalDerivative(Polynomial<T> t);
	public DifferentialForms<T> differentialForms();
	
	public Polynomial<T> resultant(Polynomial<T> t1, Polynomial<T> t2, int variable);
	
	public boolean isHomogeneous(Polynomial<T> t);
	public Polynomial<T> homogenize(Polynomial<T> t);
	public Polynomial<T> homogenize(Polynomial<T> t, int coord);
	public Polynomial<T> dehomogenize(Polynomial<T> t, int coord);
			
	public Polynomial<T> round(Polynomial<T> t, int degree);
	
	public Vector<T> asVector(Polynomial<T> t);
	public Vector<T> asVector(Polynomial<T> t, int[] degrees);
	public Polynomial<T> fromVector(Vector<T> t, int[] degrees);
}
