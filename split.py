import argparse
import csv
import random
from pathlib import Path


def split_csv(input_path, train_path, test_path, train_ratio=0.8, seed=2026):
    with input_path.open(newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    rng = random.Random(seed)
    rng.shuffle(rows)

    split_at = int(len(rows) * train_ratio)
    train_rows = rows[:split_at]
    test_rows = rows[split_at:]

    for path, output_rows in ((train_path, train_rows), (test_path, test_rows)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(output_rows)

    return len(train_rows), len(test_rows)


def main():
    parser = argparse.ArgumentParser(description="Randomly split yacht data into 80/20 train/test CSV files.")
    parser.add_argument("--input", default="data/yacht_hydro.csv", help="Input CSV path")
    parser.add_argument("--train", default="data/yacht_train.csv", help="Output training CSV path")
    parser.add_argument("--test", default="data/yacht_test.csv", help="Output test CSV path")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducible splitting")
    args = parser.parse_args()

    train_n, test_n = split_csv(
        Path(args.input),
        Path(args.train),
        Path(args.test),
        seed=args.seed,
    )
    print(f"Wrote {train_n} training rows and {test_n} test rows.")


if __name__ == "__main__":
    main()
