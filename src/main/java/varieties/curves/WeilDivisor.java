package varieties.curves;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.helper.AbstractElement;
import fields.interfaces.Element;
import varieties.projective.ProjectivePoint;

public class WeilDivisor<T extends Element<T>> extends AbstractElement<WeilDivisor<T>> {
	private Map<ProjectivePoint<T>, Integer> points;
	private int degree;

	public WeilDivisor(Map<ProjectivePoint<T>, Integer> points) {
		this.points = new TreeMap<ProjectivePoint<T>, Integer>();
		this.degree = 0;
		for (ProjectivePoint<T> p : points.keySet()) {
			if (points.get(p) != 0)
				this.points.put(p, points.get(p));
			this.degree += points.get(p);
		}
	}

	public WeilDivisor(List<ProjectivePoint<T>> zeroes, List<ProjectivePoint<T>> poles) {
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
		if (!(o instanceof WeilDivisor))
			return false;
		@SuppressWarnings("unchecked")
		WeilDivisor<T> d = (WeilDivisor<T>) o;
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
	public int compareTo(WeilDivisor<T> div) {
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