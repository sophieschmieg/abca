package fields.helper;

import java.math.BigInteger;

import fields.interfaces.Algebra;
import fields.interfaces.Element;

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
}
