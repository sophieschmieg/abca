package fields.helper;

import java.math.BigInteger;
import java.util.List;

import fields.interfaces.Algebra;
import fields.interfaces.Element;
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
	public List<T> nonTrivialCombination(List<S> s) {
		return nonTrivialCombinations(s).get(0);
	}

}
