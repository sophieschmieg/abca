package fields.interfaces;

import java.math.BigInteger;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;

import fields.exceptions.InfinityException;
import fields.polynomials.GenericUnivariatePolynomial;
import fields.polynomials.Monomial;
import fields.vectors.Vector;

public interface UnivariatePolynomialRing<T extends Element<T>> extends PolynomialRing<T> {

	Exactness exactness();

	Ring<T> getRing();

	int numberOfVariables();

	Comparator<Monomial> getComparator();

	BigInteger characteristic();

	boolean isFinite();

	BigInteger getNumberOfElements() throws InfinityException;

	boolean isIntegral();

	FactorizationResult<Polynomial<T>> uniqueFactorization(Polynomial<T> t);

	boolean isFree();

	UnivariatePolynomial<T> zero();

	UnivariatePolynomial<T> one();

	UnivariatePolynomial<T> getVar();

	Polynomial<T> getVar(int i);

	UnivariatePolynomial<T> getVarPower(int power);

	UnivariatePolynomial<T> getVarPower(int i, int power);

	UnivariatePolynomial<T> getEmbedding(T t);

	UnivariatePolynomial<T> getPolynomial(List<T> coefficients);

	UnivariatePolynomial<T> getPolynomial(@SuppressWarnings("unchecked") T... coefficients);

	UnivariatePolynomial<T> getPolynomial(Map<Monomial, T> t);

	UnivariatePolynomial<T> getEmbedding(T t, int exponent);

	UnivariatePolynomial<T> getEmbedding(T t, int[] exponents);

	UnivariatePolynomial<T> getLinear(List<T> coeff);

	UnivariatePolynomial<T> toUnivariate(Polynomial<T> t);

	boolean isHomogeneous(Polynomial<T> t);

	PolynomialRing<T> addVariableWithElimination(int shift);

	UnivariatePolynomial<T> getEmbedding(Polynomial<T> t, int[] map);

	<S extends Element<S>> UnivariatePolynomial<T> getEmbedding(Polynomial<S> t, MathMap<S, T> mapping);

	Polynomial<T> getRandomElement();

	UnivariatePolynomial<T> getRandomElement(int degree);

	SortedSet<Monomial> monomials(int degree);

	Monomial getMonomial(int[] exponents);

	GenericUnivariatePolynomial<T> add(Polynomial<T> s1, Polynomial<T> s2);

	UnivariatePolynomial<T> negative(Polynomial<T> s);

	UnivariatePolynomial<T> multiply(Polynomial<T> t1, Polynomial<T> t2);

	UnivariatePolynomial<T> multiplyPower(int power, Polynomial<T> t);

	UnivariatePolynomial<T> multiply(T t1, Polynomial<T> t2);

	UnivariatePolynomial<T> multiplyPower(int power, T t1, Polynomial<T> t2);

	UnivariatePolynomial<T> divideScalar(Polynomial<T> dividend, T divisor);

	boolean isZeroDivisor(Polynomial<T> t);

	boolean isUnit(Polynomial<T> t);

	UnivariatePolynomial<T> inverse(Polynomial<T> t);

	boolean isDivisible(Polynomial<T> dividend, Polynomial<T> divisor);

	GeneralQuotientAndRemainderResult<T> generalQuotientAndRemainder(Polynomial<T> dividend,
			List<Polynomial<T>> divisors);

	UnivariatePolynomial<T> pseudoRemainder(Polynomial<T> dividend, Polynomial<T> divisor);

	BigInteger euclidMeasure(Polynomial<T> t);

	Iterable<Polynomial<T>> getUnits();

	Map<Polynomial<T>, Integer> squareFreeFactorization(Polynomial<T> t);

	List<Polynomial<T>> getAlgebraGenerators();

	Vector<T> asVector(Polynomial<T> s);

	Vector<T> asVector(Polynomial<T> s, int degree);

	UnivariatePolynomial<T> fromVector(Vector<T> vector);

	boolean isLinearIndependent(List<Polynomial<T>> s);

	List<List<T>> nonTrivialCombinations(List<Polynomial<T>> s);

	Iterator<Polynomial<T>> iterator();

	Iterator<UnivariatePolynomial<T>> polynomials(int degree);

	Iterable<UnivariatePolynomial<T>> polynomialSet(int degree);

	Iterator<UnivariatePolynomial<T>> monicPolynomials(int degree);

	Iterable<UnivariatePolynomial<T>> monicPolynomialSet(int degree);

	T evaluate(Polynomial<T> t, List<T> ts);

	Polynomial<T> partiallyEvaluate(Polynomial<T> t, List<T> ts);

	UnivariatePolynomial<T> substitute(Polynomial<T> t, List<Polynomial<T>> values);

	boolean canBeNormalized(Polynomial<T> t);

	UnivariatePolynomial<T> normalize(Polynomial<T> t);

	T content(Polynomial<T> t);

	UnivariatePolynomial<T> contentFree(Polynomial<T> t);

	UnivariatePolynomial<T> derivative(Polynomial<T> t, int variable);

	UnivariatePolynomial<T> derivative(Polynomial<T> t);

	Vector<T> trace(Polynomial<T> t, int degree);

	T resultant(Polynomial<T> t1, Polynomial<T> t2);

	Polynomial<T> resultant(Polynomial<T> t1, Polynomial<T> t2, int variable);

	PolynomialRing<T> eliminateVariable();

	UnivariatePolynomial<Polynomial<T>> asUnivariatePolynomial(Polynomial<T> t, int variable);

	Polynomial<T> fromUnivariatePolynomial(UnivariatePolynomial<Polynomial<T>> t, int variable);

	UnivariatePolynomial<T> gcd(Polynomial<T> t1, Polynomial<T> t2);

	ExtendedEuclideanResult<Polynomial<T>> extendedEuclidean(Polynomial<T> t1, Polynomial<T> t2);

	T discriminant(Polynomial<T> t1);

	public static class ExtendedResultantResult<T extends Element<T>> {
		private final T resultant;
		private final UnivariatePolynomial<T> gcd;
		private final UnivariatePolynomial<T> coeff1;
		private final UnivariatePolynomial<T> coeff2;

		public ExtendedResultantResult(T resultant, UnivariatePolynomial<T> gcd, UnivariatePolynomial<T> coeff1,
				UnivariatePolynomial<T> coeff2) {
			this.resultant = resultant;
			this.gcd = gcd;
			this.coeff1 = coeff1;
			this.coeff2 = coeff2;
		}

		
		public T getResultant() {
			return resultant;
		}

		public UnivariatePolynomial<T> getGcd() {
			return gcd;
		}

		public UnivariatePolynomial<T> getCoeff1() {
			return coeff1;
		}

		public UnivariatePolynomial<T> getCoeff2() {
			return coeff2;
		}

	}

	ExtendedResultantResult<T> extendedResultant(Polynomial<T> t1, Polynomial<T> t2);
	
	UnivariatePolynomial<T> round(Polynomial<T> t, int degree);
}