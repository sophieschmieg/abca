package varieties.curves;

import java.util.List;

import varieties.ProjectivePoint;
import varieties.Variety;
import varieties.curves.DivisorGroup.Divisor;
import fields.Element;
import fields.Polynomial;
import fields.RationalFunction;

public interface SmoothCurve<T extends Element> extends Variety<T> {
	public int getEmbeddingDimension();
	public boolean isProjective();
	public List<Polynomial<T>> getCotangentSpace(ProjectivePoint<T> p);
	public List<Polynomial<T>> getTangentSpace(ProjectivePoint<T> p);
	public List<RationalFunction<T>> getRiemannRochSpace(Divisor<T> div);
	public boolean isPrincipal(Divisor<T> div);
}
