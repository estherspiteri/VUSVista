import React, { useState } from "react";
import styles from "./vus-upload-page.module.scss";
import Text from "./../../atoms/text/text";
import VusUploadField from "./vus-upload-field/vus-upload-field";
import Button from "../../atoms/button/button";
import PhenotypeSelection from "../sample-page/phenotype-selection/phenotype-selection";

type VusUploadPageProps = {};

//TODO: In the end show variant summary - let user confirm on cancel submittion to continue editing
const VusUploadPage: React.FunctionComponent<VusUploadPageProps> = (
  props: VusUploadPageProps
) => {
  const [chromosome, setChromosome] = useState<string | undefined>(undefined);
  const [chromosomePosition, setChromosomePosition] = useState<number>(0);
  const [refAllele, setRefAllele] = useState<string | undefined>(undefined);
  const [altAllele, setAltAllele] = useState<string | undefined>(undefined);
  const [genotype, setGenotype] = useState<string | undefined>(undefined);

  return (
    <div className={styles["vus-upload-page-container"]}>
      <div className={styles.title}>VUS Upload</div>
      <div className={styles.description}>
        <p>Input variant details in the form below.</p>
      </div>

      <div className={styles["fields-wrapper"]}>
        {/** Chromsome & Chromsome Position*/}
        <VusUploadField
          title="Chromosome"
          isOpenByDefault={true}
          showCheckMark={chromosome !== undefined && chromosomePosition !== undefined}
        >
          <div>
            <div>
              <span>Select chromosome: </span>
              <div className={styles.chromosomes}>
                {[
                  "1",
                  "2",
                  "3",
                  "4",
                  "5",
                  "6",
                  "7",
                  "8",
                  "9",
                  "10",
                  "11",
                  "12",
                  "13",
                  "14",
                  "15",
                  "16",
                  "17",
                  "18",
                  "19",
                  "20",
                  "21",
                  "22",
                  "X",
                  "Y",
                ].map((c) => (
                  <div
                    className={`${styles.chromosome} ${
                      c === chromosome ? styles["selected-chromosome"] : ""
                    }`}
                    onClick={() => setChromosome(c)}
                  >
                    {c}
                  </div>
                ))}
              </div>
            </div>

            <div>
              <span>Chromosome Position: </span>
              <Text
                value={chromosomePosition}
                type="number"
                min={0}
                className={styles["chromosome-position"]}
                onChange={(e) =>
                  setChromosomePosition(parseInt(e.currentTarget.value))
                }
              />
            </div>
          </div>
        </VusUploadField>

        {/** Reference & Alt Alleles */}
        {/**TODO: check on attempt to submit if allele consists of just GACT*/}
        <VusUploadField title="Alleles">
          <div className={styles["allele-wrapper"]}>
            <div>
              <span>Reference allele:</span>
              <Text
                value={refAllele}
                onChange={(e) => setRefAllele(e.currentTarget.value)}
              />
            </div>
            <div>
              <span>Alternate allele:</span>
              <Text
                value={altAllele}
                onChange={(e) => setAltAllele(e.currentTarget.value)}
              />
            </div>
          </div>
        </VusUploadField>

        {/** Genotype */}
        <VusUploadField title="Genotype">
          <div className={styles.genotypes}>
            {["heterozygous", "homozygous"].map((g) => {
              return (
                <div
                  className={`${styles.genotype} ${
                    genotype === g ? styles["selected-genotype"] : ""
                  }`}
                  onClick={() => setGenotype(g)}
                >
                  {g}
                </div>
              );
            })}
          </div>
        </VusUploadField>

        {/** Gene */}
        <VusUploadField title="Gene"></VusUploadField>

        {/** Samples */}
        <VusUploadField title="Samples">
          Input samples Input Samples' phenotypes
          <PhenotypeSelection />
        </VusUploadField>
      </div>

      <Button text="Save variant" />
    </div>
  );
};

export default VusUploadPage;
