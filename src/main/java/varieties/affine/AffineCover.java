package varieties.affine;

import java.util.List;

import fields.interfaces.Element;

public class AffineCover<T extends Element<T>> {
	private List<AffineScheme<T>> cover;
	private List<List<AffineMorphism<T>>> atlas;
	
	public AffineCover(List<AffineScheme<T>> cover, List<List<AffineMorphism<T>>> atlas) {
		this.cover = cover;
		this.atlas = atlas;
		if (cover.size() != atlas.size()) {
			throw new ArithmeticException("Not enough intersections!");
		}
		for (int i = 0; i < cover.size(); i++) {
			if (atlas.get(i).size() != cover.size()) {
				throw new ArithmeticException("Not enough intersections!");
			}
			for (int j = 0; j < i; j++) {
				if (!atlas.get(i).get(j).getDomain().equals(atlas.get(j).get(i).getDomain())) {
					throw new ArithmeticException("Not the same intersection!");
				}
				if (!atlas.get(i).get(j).getRange().equals(cover.get(i))) {
					throw new ArithmeticException("Not the correct range!");
				}
				if (!atlas.get(j).get(i).getRange().equals(cover.get(j))) {
					throw new ArithmeticException("Not the correct range!");
				}
				if (!atlas.get(i).get(j).isInjective()) {
					throw new ArithmeticException("Not injective!");
				}
				if (!atlas.get(j).get(i).isInjective()) {
					throw new ArithmeticException("Not injective!");
				}
			}
		}
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
}
