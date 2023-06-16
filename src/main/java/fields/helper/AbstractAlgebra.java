package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import fields.interfaces.Algebra;
import fields.interfaces.Element;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

public abstract class AbstractAlgebra<T extends Element<T>, S extends Element<S>> extends AbstractRing<S>
		implements Algebra<T, S> {
	public AbstractAlgebra() {
		super(true);
	}

	public AbstractAlgebra(boolean generateUnivariatePolynomialRing) {
		super(generateUnivariatePolynomialRing);
	}

	@Override
	public S scalarMultiply(T t, S s) {
		return multiply(getEmbedding(t), s);
	}
	
	@Override
	public S scalarMultiply(int n, S s) {
		return scalarMultiply(getRing().getInteger(n), s);
	}

	@Override
	public S scalarMultiply(BigInteger n, S s) {
		return scalarMultiply(getRing().getInteger(n), s);
	}

	@Override
	public S scalarMultiply(T t1, T t2, S s) {
		return scalarMultiply(getRing().multiply(t1, t2), s);
	}

	@Override
	public S scalarMultiply(int n, T t, S s) {
		return scalarMultiply(getRing().multiply(n, t), s);
	}

	@Override
	public S scalarMultiply(BigInteger n, T t, S s) {
		return scalarMultiply(getRing().multiply(n, t), s);
	}

	@Override
	public S scalarMultiply(int n, T t1, T t2, S s) {
		return scalarMultiply(getRing().multiply(n, t1, t2), s);
	}

	@Override
	public S scalarMultiply(BigInteger n, T t1, T t2, S s) {
		return scalarMultiply(getRing().multiply(n, t1, t2), s);
	}

	@Override
	public S fromVector(Vector<T> vector) {
		S result = zero();
		List<S> generators = getModuleGenerators();
		for (int i = 0; i < generators.size(); i++) {
			result = add(result, scalarMultiply(vector.get(i+1), generators.get(i)));
		}
		return result;
	}

	@Override
	public boolean isTorsionFree() {
		return annihilator().equals(getRing().getZeroIdeal());
	}
	
	private Matrix<T> asMatrix(List<S> m) {
		List<Vector<T>> result = new ArrayList<>();
		for (S s : m) {
			result.add(asVector(s));
		}
		return Matrix.fromColumns(result);
	}
	
	@Override
	public boolean isSubModuleMemberModule(List<S> m, S b) {
		Matrix<T> matrix = asMatrix(m);
		return isSubModuleMemberModule(matrix.getModule(getRing()), matrix, asVector(b));
	}
	
	@Override
	public boolean isSubModuleMemberModule(MatrixModule<T> module, Matrix<T> m, Vector<T> b) {
		return getRing().isSubModuleMember(module, m, b);
	}
	
	@Override
	public S asSubModuleMemberModule(List<S> m, S b) {
		Matrix<T> matrix = asMatrix(m);
		return fromVector(asSubModuleMemberModule(matrix.getModule(getRing()), matrix, asVector(b)));
	}
	
	@Override
	public Vector<T> asSubModuleMemberModule(MatrixModule<T> module, Matrix<T> m, Vector<T> b) {
		return getRing().asSubModuleMember(module, m, b);
	}
	
	@Override
	public List<Vector<T>> syzygyProblemModule(List<S> m) {
		Matrix<T> matrix = asMatrix(m);
		return syzygyProblemModule(matrix.getModule(getRing()), matrix);
	}
	
	@Override
	public List<Vector<T>> syzygyProblemModule(MatrixModule<T> module, Matrix<T> m) {
		return getRing().syzygyProblem(module, m);
	}
	
	@Override
	public List<S> simplifySubModuleGeneratorsModule(List<S> m) {
		Matrix<T> matrix = asMatrix(m);
		List<Vector<T>> simplified = simplifySubModuleGeneratorsModule(matrix.getModule(getRing()), matrix);
		List<S> result = new ArrayList<>();
		for (Vector<T> vector : simplified) {
			result.add(fromVector(vector));
		}
		return result;
	}
	
	@Override
	public List<Vector<T>> simplifySubModuleGeneratorsModule(MatrixModule<T> module, Matrix<T> m) {
		return getRing().simplifySubModuleGenerators(module, m);
	}

}
