package fields.interfaces;

import java.util.List;

import fields.helper.FieldAutomorphism;
import fields.helper.FieldEmbedding;
import fields.helper.GaloisGroup;
import fields.vectors.FiniteVectorSpace;
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
public interface FieldExtension<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, Ext extends FieldExtension<T, S, Ext>>
		extends Field<S>, VectorSpace<T, S>, AlgebraicRingExtension<T, S, Ext> {
	Field<T> getBaseField();

	int degree();

	int degreeOver(FieldEmbedding<T, S, Ext> base);

	S getEmbedding(T t);

	MathMap<T, S> getEmbeddingMap();

	/**
	 * alpha() returns an element such that L = K[alpha].
	 */
	S alpha();

	UnivariatePolynomial<T> minimalPolynomial();

	UnivariatePolynomial<S> minimalPolynomialOver(FieldEmbedding<T, S, Ext> base);

	GaloisGroup<T, S, Ext> galoisGroup();

	boolean isSeparable();

	boolean isNormal();

	boolean isFiniteExtension();

	boolean isGalois();

	boolean isCyclic();

	FieldAutomorphism<T, S, Ext> cyclicGenerator();

	/**
	 * Only allowed if totalDegree() == 1
	 */
	T asBaseFieldElement(S s);

	MathMap<S, T> asBaseFieldElementMap();

	FiniteVectorSpace<T> asVectorSpace();

	Vector<T> asVector(S s);

	MathMap<S, Vector<T>> asVectorMap();

	Vector<S> asVectorOver(S s, FieldEmbedding<T, S, Ext> base);

	MatrixAlgebra<T> matrixAlgebra();

	Matrix<T> asMatrix(S s);

	Matrix<S> asMatrixOver(S s, FieldEmbedding<T, S, Ext> base);

	UnivariatePolynomial<T> asPolynomial(S s);

	S fromPolynomial(UnivariatePolynomial<T> s);

	PolynomialRing<T> genericPolynomialRing();

	FreeModule<Polynomial<T>> genericFreeModule();

	Vector<Polynomial<T>> genericVector();

	Matrix<Polynomial<T>> genericMatrix();

	Polynomial<T> genericNorm();

	List<S> conjugates(S s);

	S hilbert90(S s);

	T norm(S s);

	T trace(S s);

	S normOver(S s, FieldEmbedding<T, S, Ext> base);

	S traceOver(S s, FieldEmbedding<T, S, Ext> base);
	
	T traceForm(S s1, S s2);
	
	S traceFormOver(S s1, S s2, FieldEmbedding<T, S, Ext> base);
	
	BilinearMap<S, T> traceForm();
	
	BilinearMap<S, S> traceFormOver(FieldEmbedding<T, S, Ext> base);
	
	Matrix<T> traceFormMatrix();
	
	Matrix<S> traceFormMatrixOver(FieldEmbedding<T, S, Ext> base);

	UnivariatePolynomial<T> minimalPolynomial(S s);

	UnivariatePolynomial<S> minimalPolynomialOver(S s, FieldEmbedding<T, S, Ext> base);

	Extension<S, T, S, Ext> getExtension(UnivariatePolynomial<S> minimalPolynomial);

	FieldEmbedding<T, S, Ext> getEmbeddedExtension(UnivariatePolynomial<S> minimalPolynomial);

	Ext makeExtension(UnivariatePolynomial<T> minimalPolynomial);
}
