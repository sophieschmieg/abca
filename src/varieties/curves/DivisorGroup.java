package varieties.curves;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.Element;
import fields.Group;
import fields.InfinityException;
import varieties.ProjectivePoint;
import varieties.curves.DivisorGroup.Divisor;

public class DivisorGroup<T extends Element> implements Group<Divisor<T>> {

	@Override
	public Divisor<T> getRandomElement() {
		return neutral(); // Chosen by fair dice roll
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<Divisor<T>> iterator() {
		throw new InfinityException();
	}

	@Override
	public Divisor<T> neutral() {
		Map<ProjectivePoint<T>, Integer> map = Collections.emptyMap();
		return new Divisor<T>(map);
	}

	@Override
	public Divisor<T> inverse(Divisor<T> t) {
		Map<ProjectivePoint<T>, Integer> map = new TreeMap<ProjectivePoint<T>, Integer>();
		for (ProjectivePoint<T> p : t.getSupport())
			map.put(p, -t.get(p));
		return new Divisor<T>(map);
	}

	@Override
	public Divisor<T> operate(Divisor<T> t1, Divisor<T> t2) {
		Map<ProjectivePoint<T>, Integer> map = new TreeMap<ProjectivePoint<T>, Integer>();
		Set<ProjectivePoint<T>> keyset = new TreeSet<ProjectivePoint<T>>();
		keyset.addAll(t1.getSupport());
		keyset.addAll(t2.getSupport());
		for (ProjectivePoint<T> p : keyset) {
			map.put(p, t1.get(p) + t2.get(p));
		}
		return new Divisor<T>(map);
	}

	public static class Divisor<T extends Element> implements Element {
		private Map<ProjectivePoint<T>, Integer> points;
		private int degree;

		public Divisor(Map<ProjectivePoint<T>, Integer> points) {
			this.points = new TreeMap<ProjectivePoint<T>, Integer>();
			this.degree = 0;
			for (ProjectivePoint<T> p : points.keySet()) {
				if (points.get(p) != 0)
					this.points.put(p, points.get(p));
				this.degree += points.get(p);
			}
		}

		public Divisor(List<ProjectivePoint<T>> zeroes, List<ProjectivePoint<T>> poles) {
			this.points = new TreeMap<ProjectivePoint<T>, Integer>();
			for (ProjectivePoint<T> zero : zeroes) {
				if (!this.points.containsKey(zero))
					this.points.put(zero, 1);
				else
					this.points.put(zero, this.points.get(zero) + 1);
			}
			for (ProjectivePoint<T> pole : poles) {
				if (!this.points.containsKey(pole))
					this.points.put(pole, -1);
				else if (this.points.get(pole) == 1)
					this.points.remove(pole);
				else
					this.points.put(pole, this.points.get(pole) - 1);
			}
			this.degree = zeroes.size() - poles.size();
		}

		public int getDegree() {
			return degree;
		}

		public Set<ProjectivePoint<T>> getSupport() {
			return Collections.unmodifiableSet(points.keySet());
		}

		public Integer get(ProjectivePoint<T> point) {
			if (this.points.containsKey(point))
				return this.points.get(point);
			return 0;
		}

		public List<ProjectivePoint<T>> getZeroes() {
			List<ProjectivePoint<T>> zeroes = new ArrayList<ProjectivePoint<T>>();
			for (ProjectivePoint<T> p : this.points.keySet()) {
				int num = this.points.get(p);
				if (num > 0)
					for (int i = 0; i < num; i++)
						zeroes.add(p);
			}
			return Collections.unmodifiableList(zeroes);
		}

		public List<ProjectivePoint<T>> getPoles() {
			List<ProjectivePoint<T>> poles = new ArrayList<ProjectivePoint<T>>();
			for (ProjectivePoint<T> p : this.points.keySet()) {
				int num = this.points.get(p);
				if (num < 0)
					for (int i = 0; i < -num; i++)
						poles.add(p);
			}
			return Collections.unmodifiableList(poles);
		}

		public boolean isEffective() {
			for (ProjectivePoint<T> p : this.points.keySet()) {
				if (this.points.get(p) < 0)
					return false;
			}
			return true;
		}

		@Override
		public boolean equals(Object o) {
			if (!(o instanceof Divisor))
				return false;
			@SuppressWarnings("unchecked")
			Divisor<T> d = (Divisor<T>) o;
			return this.points.equals(d.points);
		}

		@Override
		public String toString() {
			StringBuffer buf = new StringBuffer();
			boolean first = true;
			for (ProjectivePoint<T> p : this.points.keySet()) {
				if (this.points.get(p) > 0) {
					if (first)
						first = false;
					else
						buf.append(" + ");
					if (this.points.get(p) > 1)
						buf.append(this.points.get(p));
					buf.append(p);
				}
			}
			for (ProjectivePoint<T> p : this.points.keySet()) {
				if (this.points.get(p) < 0) {
					buf.append(" - ");
					if (this.points.get(p) < -1)
						buf.append((-1) * this.points.get(p));
					buf.append(p);
				}
			}
			return buf.toString();
		}

		@Override
		public int compareTo(Element o) {
			@SuppressWarnings("unchecked")
			Divisor<T> div = (Divisor<T>) o;
			Set<ProjectivePoint<T>> points = new TreeSet<ProjectivePoint<T>>();
			points.addAll(this.points.keySet());
			points.addAll(div.points.keySet());
			for (ProjectivePoint<T> p : points) {
				int cmp = this.get(p).compareTo(div.get(p));
				if (cmp != 0)
					return cmp;
			}
			return 0;
		}
	}
}
