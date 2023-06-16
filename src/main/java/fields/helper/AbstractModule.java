package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.interfaces.Element;
import fields.interfaces.Group;
import fields.interfaces.Module;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

public abstract class AbstractModule<T extends Element<T>, S extends Element<S>> implements Module<T, S> {

	@Override
	public Group<S> getAdditiveGroup() {
		return new Group<S>() {

			@Override
			public Exactness exactness() {
				return AbstractModule.this.exactness();
			}

			@Override
			public S getRandomElement() {
				return AbstractModule.this.getRandomElement();
			}

			@Override
			public boolean isFinite() {
				return AbstractModule.this.isFinite();
			}

			@Override
			public BigInteger getNumberOfElements() throws InfinityException {
				return AbstractModule.this.getNumberOfElements();
			}

			@Override
			public Iterator<S> iterator() {
				return AbstractModule.this.iterator();
			}

			@Override
			public S neutral() {
				return AbstractModule.this.zero();
			}

			@Override
			public S inverse(S s) {
				return AbstractModule.this.negative(s);
			}

			@Override
			public S operate(S s1, S s2) {
				return AbstractModule.this.add(s1, s2);
			}

		};
	}
	
	@Override
	public S add(S s1, S s2, S s3) {
		return add(add(s1, s2), s3);
	}

	@Override
	public S subtract(S minuend, S subtrahend) {
		return add(minuend, negative(subtrahend));
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
