package util;

import java.util.ArrayList;
import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.MathMap;
import fields.vectors.Vector;

public class FlattenProductVectorMap<T extends Element<T>, S extends Element<S>> implements MathMap<Vector<T>, Vector<S>> {
	private MathMap<T, Vector<S>> map;
	public FlattenProductVectorMap(MathMap<T, Vector<S>> map) {
		this.map = map;
	}
		
		@Override
		public Vector<S> evaluate(Vector<T> t) {
			List<S> results = new ArrayList<>();
			for (T e : t.asList()) {
				results.addAll(map.evaluate(e).asList());
			}
			return new Vector<>(results);
		}

	}
