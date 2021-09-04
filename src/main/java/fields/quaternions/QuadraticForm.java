package fields.quaternions;

import java.util.ArrayList;
import java.util.List;

import fields.interfaces.BilinearMap;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.VectorSpace;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.SubVectorSpace;
import fields.vectors.Vector;

public class QuadraticForm<T extends Element<T>, S extends Element<S>> implements MathMap<S, T>{
	private PolynomialRing<T> polynomialRing;
	private Polynomial<T> asPolynomial;
	private BilinearMap<S, T> asBilinearMap;
	private Matrix<T> asMatrix;
	private VectorSpace<T, S> space;
	private Field<T> field;

	public QuadraticForm(VectorSpace<T, S> space, Polynomial<T> asPolynomial) {
		if (asPolynomial.degree() != 2) {
			throw new ArithmeticException("Not a quadratic form!");
		}
		this.asPolynomial = asPolynomial;
		this.polynomialRing = asPolynomial.getPolynomialRing();
		if (polynomialRing.getRing().equals(space.getField())) {
			throw new ArithmeticException("Fields mismatched!");
		}
		if (space.dimension() != polynomialRing.numberOfVariables()) {
			throw new ArithmeticException("Number of variables mismatched!");
		}
		this.space = space;
		this.field = space.getField();
		this.asBilinearMap = new BilinearMap<>() {

			@Override
			public T evaluate(S t1, S t2) {
				T xplusy = polynomialRing.evaluate(asPolynomial, space.asVector(space.add(t1, t2)));
				T x = polynomialRing.evaluate(asPolynomial, space.asVector(t1));
				T y = polynomialRing.evaluate(asPolynomial, space.asVector(t2));
				return field.add(xplusy, field.negative(field.add(x, y)));
			}
		};
		this.asMatrix = Matrix.fromBilinearMap(space, asBilinearMap);
	}

	public SubVectorSpace<T, S> orthogonalSpace(SubVectorSpace<T, S> subspace) {
		List<Vector<T>> asRows = new ArrayList<>();
		for (S s : subspace.getBasis()) {
			asRows.add(space.asVector(s));
		}
		Matrix<T> rowMatrix = Matrix.fromRows(asRows);
		MatrixModule<T> module = new MatrixModule<>(field, subspace.dimension(), space.dimension());
		Matrix<T> curryMatrix = module.multiply(rowMatrix, asMatrix);
		List<Vector<T>> basis = module.kernelBasis(curryMatrix);
		List<S> fromVectors = new ArrayList<>();
		for (Vector<T> b : basis) {
			fromVectors.add(space.fromVector(b));
		}
		return new SubVectorSpace<>(space, fromVectors);
	}
	
	public T discriminant() {
		return field.multiply(field.power(field.getInteger(2), -space.dimension()), space.matrixAlgebra().determinant(asMatrix));
	}
	
	@Override
	public T evaluate(S s) {
		return polynomialRing.evaluate(asPolynomial, space.asVector(s));
	}
	
	public Polynomial<T> asPolynomial() {
		return asPolynomial;
	}
}
