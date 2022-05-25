package fields.interfaces;

import java.util.Comparator;
import java.util.List;

import fields.polynomials.CoordinateRing;
import fields.polynomials.Monomial;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
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

	String getVariableName();
	
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

	CoordinateRing<T> asCoordinateRing();

	public static class PolynomialRingAsCoordinateRing<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>> {
		private PolynomialRing<S> polynomialRing;
		private CoordinateRing<T> coordinateRing;
		private MathMap<Polynomial<S>, CoordinateRingElement<T>> isomorphism;
		private MathMap<CoordinateRingElement<T>, Polynomial<S>> inverseIsomorphism;

		public PolynomialRingAsCoordinateRing(PolynomialRing<S> polynomialRing, CoordinateRing<T> coordinateRing,
				MathMap<Polynomial<S>, CoordinateRingElement<T>> isomorphism,
				MathMap<CoordinateRingElement<T>, Polynomial<S>> inverseIsomorphism) {
			this.polynomialRing = polynomialRing;
			this.coordinateRing = coordinateRing;
			this.isomorphism = isomorphism;
			this.inverseIsomorphism = inverseIsomorphism;
		}

		public PolynomialRing<S> getPolynomialRing() {
			return polynomialRing;
		}

		public CoordinateRing<T> getCoordinateRing() {
			return coordinateRing;
		}

		public MathMap<Polynomial<S>, CoordinateRingElement<T>> getIsomorphism() {
			return isomorphism;
		}

		public MathMap<CoordinateRingElement<T>, Polynomial<S>> getInverseIsomorphism() {
			return inverseIsomorphism;
		}
	}

	PolynomialRingAsCoordinateRing<T, S> asCoordinateRing(PolynomialRing<S> polynomialRing);

	PolynomialRingAsCoordinateRing<T, S> asCoordinateRing(int numberOfVariables);

	PolynomialRingAsCoordinateRing<T, S> asCoordinateRing(int numberOfVariables, Comparator<Monomial> comparator);


	public static class ExtensionCoordinateRing<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>> {
		private CoordinateRing<S> extensionCoordinateRing;
		private CoordinateRing<T> baseCoordinateRing;
		private MathMap<CoordinateRingElement<S>, CoordinateRingElement<T>> isomorphism;
		private MathMap<CoordinateRingElement<T>, CoordinateRingElement<S>> inverseIsomorphism;

		public ExtensionCoordinateRing(CoordinateRing<S> extensionCoordinateRing, CoordinateRing<T> baseCoordinateRing,
				MathMap<CoordinateRingElement<S>, CoordinateRingElement<T>> isomorphism,
				MathMap<CoordinateRingElement<T>, CoordinateRingElement<S>> inverseIsomorphism) {
			this.extensionCoordinateRing = extensionCoordinateRing;
			this.baseCoordinateRing = baseCoordinateRing;
			this.isomorphism = isomorphism;
			this.inverseIsomorphism = inverseIsomorphism;
		}

		public CoordinateRing<S> getExtensionCoordinateRing() {
			return extensionCoordinateRing;
		}

		public CoordinateRing<T> getBaseCoordinateRing() {
			return baseCoordinateRing;
		}

		public MathMap<CoordinateRingElement<S>, CoordinateRingElement<T>> getIsomorphism() {
			return isomorphism;
		}

		public MathMap<CoordinateRingElement<T>, CoordinateRingElement<S>> getInverseIsomorphism() {
			return inverseIsomorphism;
		}
	}

	ExtensionCoordinateRing<T, S> asCoordinateRing(CoordinateRing<S> coordinateRing);
}
