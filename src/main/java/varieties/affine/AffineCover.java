package varieties.affine;

import java.util.ArrayList;
import java.util.List;

import fields.interfaces.Element;
import varieties.affine.AffineScheme.AffineSimplification;

public class AffineCover<T extends Element<T>> {
	private List<AffineScheme<T>> cover;
	private List<List<AffineMorphism<T>>> atlas;

	public AffineCover(List<AffineScheme<T>> cover, List<List<AffineMorphism<T>>> atlas) {
		this(cover, atlas, true);
	}

	public AffineCover(List<AffineScheme<T>> cover, List<List<AffineMorphism<T>>> atlas, boolean checkAtlas) {
		this.cover = cover;
		this.atlas = atlas;
		if (cover.size() != atlas.size()) {
			throw new ArithmeticException("Not enough intersections!");
		}
		for (int i = 0; i < cover.size(); i++) {
			if (atlas.get(i).size() != cover.size()) {
				throw new ArithmeticException("Not enough intersections!");
			}
			if (checkAtlas) {
				for (int j = 0; j <= i; j++) {
					if (!atlas.get(i).get(j).getDomain().equals(atlas.get(j).get(i).getDomain())) {
						throw new ArithmeticException("Not the same intersection!");
					}
					if (!atlas.get(i).get(j).getRange().equals(cover.get(i))) {
						throw new ArithmeticException("Not the correct range!");
					}
					if (!atlas.get(j).get(i).getRange().equals(cover.get(j))) {
						throw new ArithmeticException("Not the correct range!");
					}
					if (!atlas.get(i).get(j).isOpenImmersion()) {
						throw new ArithmeticException("Not an open immersion: i = " + i + ", j = " + j + "!");
					}
					if (!atlas.get(j).get(i).isOpenImmersion()) {
						throw new ArithmeticException("Not an open immersion: j = " + j + ", i = " + i + "!");
					}
				}
			}
		}
	}

	public int size() {
		return cover.size();
	}

	public AffineScheme<T> get(int i) {
		return cover.get(i);
	}

	public List<AffineScheme<T>> getCover() {
		return cover;
	}

	public AffineScheme<T> getIntersection(int i, int j) {
		return atlas.get(i).get(j).getDomain();
	}

	public AffineMorphism<T> getIntersectionEmbedding(int i, int other) {
		return atlas.get(i).get(other);
	}

	public AffineCover<T> simplifyIntersections() {
		List<AffineScheme<T>> simplerCover = new ArrayList<>();
		simplerCover.addAll(cover);
		List<List<AffineMorphism<T>>> simplerAtlas = new ArrayList<>();
		for (int i = 0; i < cover.size(); i++) {
			simplerAtlas.add(new ArrayList<>());
			for (int j = 0; j < cover.size(); j++) {
				simplerAtlas.get(i).add(atlas.get(i).get(j));
			}
		}
		for (int i = 0; i < cover.size(); i++) {
			for (int j = 0; j <= i; j++) {
				AffineSimplification<T> simplification = atlas.get(i).get(j).getDomain().simplify();
				simplerAtlas.get(i).set(j,
						AffineMorphism.composition(simplification.getInverse(), atlas.get(i).get(j)));
				simplerAtlas.get(j).set(i,
						AffineMorphism.composition(simplification.getInverse(), atlas.get(j).get(i)));
			}
		}
		return new AffineCover<>(simplerCover, simplerAtlas);
	}

	public boolean equals(Object o) {
		if (!(o instanceof AffineCover<?>)) {
			return false;
		}
		@SuppressWarnings("unchecked")
		AffineCover<T> other = (AffineCover<T>) o;
		return atlas.equals(other.atlas);
	}
}
