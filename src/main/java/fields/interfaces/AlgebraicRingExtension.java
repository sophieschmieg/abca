package fields.interfaces;

import java.util.List;

import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;

/**
 * An interface to describe an algebraic field extension L / K / P with P being
 * a common smallest subfield.
 *
 * @param <T>
 * @param <S>
 */
public interface AlgebraicRingExtension<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, Ext extends AlgebraicRingExtension<T, S, Ext>>
		extends Algebra<T, S> {
	int degree();

	S getEmbedding(T t);

	MathMap<T, S> getEmbeddingMap();

	/**
	 * alpha() returns an element such that L = K[alpha].
	 */
	S alpha();

	List<Ext> asIrreducibleProduct();

	List<S> asIrreducibleProductElement(S s);

	S fromIrreducibleProductElement(List<S> s);

	UnivariatePolynomial<T> minimalPolynomial();

	/**
	 * Only allowed if totalDegree() == 1
	 */
	T asBaseFieldElement(S s);

	MathMap<S, T> asBaseFieldElementMap();

	FreeModule<T> asFreeModule();

	Vector<T> asVector(S s);

	MathMap<S, Vector<T>> asVectorMap();

	MatrixAlgebra<T> matrixAlgebra();

	Matrix<T> asMatrix(S s);

	UnivariatePolynomial<T> asPolynomial(S s);

	S fromPolynomial(UnivariatePolynomial<T> s);

	PolynomialRing<T> genericPolynomialRing();

	FreeModule<Polynomial<T>> genericFreeModule();

	Vector<Polynomial<T>> genericVector();

	Matrix<Polynomial<T>> genericMatrix();

	UnivariatePolynomial<T> minimalPolynomial(S s);
}
