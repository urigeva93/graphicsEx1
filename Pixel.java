
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class Pixel {
	int x;
	int y;

	public Pixel(int x, int y) {
		this.x = x;
		this.y = y;
	}

	public Pixel(Pixel p) {
		this.x = p.x;
		this.y = p.y;
	}

	public int getX() {
		return x;
	}

	public void setX(int x) {
		this.x = x;
	}

	public int getY() {
		return y;
	}

	public void setY(int y) {
		this.y = y;
	}

	private static List<Pixel> removeUnrealPixels(List<Pixel> pixelsList, int w, int h) {
		Iterator<Pixel> it = pixelsList.iterator();

		List<Pixel> pixelsToRemove = new ArrayList<Pixel>();

		// First collect all nodes to remove
		while (it.hasNext()) {
			Pixel p = it.next();
			if ((p.x < 0) || (p.y < 0) || (p.x > h - 1) || (p.y > w - 1)) {
				pixelsToRemove.add(p);
			}
		}

		// Then Remove it from neighborsList
		it = pixelsToRemove.iterator();
		while (it.hasNext()) {
			Pixel p = it.next();
			pixelsList.remove(p);
		}
		return pixelsList;

	}

	public List<Pixel> getNeighbors(int w, int h)

	{
		List<Pixel> neighborsList = new ArrayList<Pixel>();

		neighborsList.add(new Pixel(this.x - 1, this.y - 1));
		neighborsList.add(new Pixel(this.x - 1, this.y));
		neighborsList.add(new Pixel(this.x - 1, this.y + 1));
		neighborsList.add(new Pixel(this.x, this.y - 1));
		neighborsList.add(new Pixel(this.x, this.y + 1));
		neighborsList.add(new Pixel(this.x + 1, this.y - 1));
		neighborsList.add(new Pixel(this.x + 1, this.y));
		neighborsList.add(new Pixel(this.x + 1, this.y + 1));

		neighborsList = removeUnrealPixels(neighborsList, w, h);
		return neighborsList;
	}

	public List<Pixel> getEnthropyMembers(int w, int h) {
		List<Pixel> neighborsList = new ArrayList<Pixel>();
		for (int i = this.x - 4; i <= this.x + 4; i++) {
			for (int j = this.y - 4; j <= this.y + 4; j++) {
				Pixel p = new Pixel(i, j);
				neighborsList.add(p);
			}
		}
		neighborsList = removeUnrealPixels(neighborsList, w, h);
		return neighborsList;
	}
}
