# crs-path-inspection-recount

論文「OWASP CRSのパス検査拡張に伴う追加遮断とルール単位の原因特定」(空口・鈴木・岡崎)の再集計資料です。
匿名化した観測記録と、表3〜7の値を数え直すスクリプトを置いています。

Anonymized observations and a recount script for Tables 3–7 of our paper on URI-path inspection in OWASP CRS 4.30.0.

## 使い方

```
python recount.py
```

Python 3.9 以降で動きます。追加のパッケージは要りません。表ごとの値を JSON で出力し、`expected.json`(論文に載せた値)と一致するかを確かめます。HTTP 要求は送りません。

## data/

| ファイル | 内容 | 論文の対応箇所 |
|---|---|---|
| `main_sample.csv` | 主標本 15,561 パスの判定(PL1〜4、G0/G9/G31、PL3・4 は G31P も) | 表3、表7(a) |
| `removal_checks.csv` | PL3 で G31 が追加遮断した 41 パス(主標本 23、補助標本 18)について、パス検査を外したときの判定 | 表4 |
| `sites.csv` | 用途確認サイト 30 件のパスの判定 | 4.1.5項、表7(b) |
| `regression_inputs.csv` | 回帰試験由来の選別群・無作為群(各 24 入力)の予測と判定 | 表5 |
| `detailed_inputs.csv` | 詳細確認の 4 入力の判定(304 送信) | 表6 |
| `random_group_candidates.csv` | 無作為群の候補 657 件と、除外・抽出の結果 | 3.2.3項 |

判定の列は HTTP 応答コードです(200 が通過、403 が遮断)。`main_sample.csv` は 1 行がパスと PL の組、`regression_inputs.csv` は 1 行が入力・PL・構成の組です。

- `G0` `G9` `G31`:932 群のうちパスを検査するルールが 0 本、9 本、31 本の構成。`G31P` は G31 の 932237 から who と es を除いた構成。
- `G31_minus_932237` など:G31 からそのルールのパス検査だけを外した構成。
- `cluster`:抽出元の Common Crawl 索引ブロック(番号は付け替えてあります)。ブロック単位の再標本化に使います。
- `path_rules_*`:パスに一致した 932 群のルール番号。
- `G31P` が空欄の行は測定していない条件です(PL1・2)。
- `sites.csv` の `path_string` は、同じパス文字列のサイトをまとめるための番号です。パスそのものではありません。
- `random_group_candidates.csv` の `source_rule` と `test_id` は CRS 4.30.0 の回帰試験の番号です。`path_key` は、選別群と同じパスになる候補を除くための照合キー(符号化後のパスの SHA-256 の先頭 16 桁)で、`regression_inputs.csv` にも同じ形で入れてあります。

`main_sample.csv` の行の順序は変えないでください。再標本化の結果が行の順序に依存します。

## 含めていないもの

ホスト名、URL、パスの文字列、ページ本文、WAF に送った要求は含めていません。そのため、この資料でできるのは保存した判定の再集計までで、WAF への再送信はできません。

## ファイルの確認

```
sha256sum -c SHA256SUMS
```

## ライセンス

スクリプトは MIT License、データは CC BY 4.0 です。
