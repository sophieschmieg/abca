package util;

import fields.interfaces.Element;
import fields.interfaces.MathMap;

public class ConCatMap<T extends Element<T>, S extends Element<S>, U extends Element<U>> implements MathMap<T, U> {
	private MathMap<T, S> first;
	private MathMap<S, U> second;

	public ConCatMap(MathMap<T, S> first, MathMap<S, U> second) {
		this.first = first;
		this.second = second;
	}

	@Override
	public U evaluate(T t) {
		return second.evaluate(first.evaluate(t));
	}

}
