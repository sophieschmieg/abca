package varieties;

import fields.Element;
import fields.Field;
import fields.MathSet;

public interface Variety<T extends Element> extends MathSet<ProjectivePoint<T>>{
	public Field<T> getField();
	public boolean hasRationalPoint(ProjectivePoint<T> p);
}
