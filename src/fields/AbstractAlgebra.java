package fields;

import java.math.BigInteger;
import java.util.List;

public abstract class AbstractAlgebra<T extends Element, S extends Element> extends AbstractRing<S>
		implements Module<T, S> {
	public abstract S getEmbedding(T t);
	public abstract List<S> getGenerators();
	
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
