import React, { useState } from "react";
import styles from "./vus-info.module.scss";
import { IVus } from "../../../models/view-vus.model";
import VariantSummary from "../../shared/variant-summary/variant-summary";
import Icon from "../../../atoms/icons/icon";
import { openInNewWindow } from "../../../helpers/open-links";
import Button from "../../../atoms/button/button";
import { IAcmgRule } from "../../../models/acmg-rule.model";
import AcmgRulesEdit from "../../sample-page/acmg-rules-edit/acmg-rules-edit";
import { VusService } from "../../../services/vus/vus.service";
import AcmgRuleInfo from "../../sample-page/acmg-rule-info/acmg-rule-info";

type VusInfoProps = {
  vus: IVus;
  acmgRules: IAcmgRule[];
  vusService?: VusService;
};

const VusInfo: React.FunctionComponent<VusInfoProps> = (
  props: VusInfoProps
) => {
  const [acmgRuleHover, setAcmgRuleHover] = useState<number | undefined>(
    undefined
  );
  const [isAcmgEditMenuOpen, setIsAcmgEditMenuOpen] = useState(false);

  return (
    <div className={styles["vus-info-container"]}>
      <div className={styles["vus-info"]}>
        <div className={styles["vus-summary"]}>
          <VariantSummary
            variant={{
              id: props.vus.id,
              chromosome: props.vus.chromosome,
              chromosomePosition: props.vus.chromosomePosition,
              gene: props.vus.gene,
              refAllele: props.vus.refAllele,
              altAllele: props.vus.altAllele,
            }}
          />
          {props.vus.rsidDbsnpVerified && (
            <Button
              text="View Publications"
              icon="publication"
              onClick={() =>
                openInNewWindow(`/publication-view/${props.vus.id}`)
              }
            />
          )}
        </div>
        <div className={styles["top-container"]}>
          {/** General information */}
          <div className={styles.information}>
            <div className={styles.info}>
              <div className={styles["info-title"]}>Type:</div>
              {props.vus.type}
            </div>

            <div className={styles.info}>
              <div className={styles["info-title"]}>Homozygotes:</div>
              {props.vus.numHomozygous}
            </div>

            <div className={styles.info}>
              <div className={styles["info-title"]}>Heterozygotes:</div>
              {props.vus.numHeterozygous}
            </div>
          </div>
        </div>

        {/*TODO: Next to is RSID verified do info icon - on hover show what info was compared. Same for clinvar*/}
        {/** External References */}
        <div className={styles["external-ref"]}>
          {/** Clinvar */}
          <div className={styles["info-container"]}>
            <div className={styles["external-ref-title-container"]}>
              <p className={styles["info-title"]}>
                {`ClinVar${
                  props.vus.clinvarClassification?.length > 0 &&
                  !props.vus.rsidDbsnpVerified
                    ? " of suggested dbSNP RSID"
                    : ""
                }:`}
              </p>
              {(props.vus.clinvarClassification?.length > 0 ||
                props.vus.clinvarErrorMsg?.length > 0) && (
                <Icon
                  name="external-link"
                  className={`${styles["external-link"]} ${styles.clinvar} ${
                    props.vus.clinvarClassification?.length === 0
                      ? styles.disabled
                      : props.vus.rsid?.length > 0 &&
                        !props.vus.rsidDbsnpVerified
                      ? styles["unverified-rsid"]
                      : styles.active
                  }`}
                  onClick={(e) => {
                    if (props.vus.clinvarErrorMsg?.length === 0) {
                      e.stopPropagation();
                      openInNewWindow(
                        `https://www.ncbi.nlm.nih.gov/clinvar/variation/${props.vus.clinvarUid}`
                      );
                    }
                  }}
                />
              )}
            </div>
            <div
              className={`${styles.info} ${
                (!props.vus.clinvarClassification ||
                  props.vus.clinvarClassification?.length === 0) &&
                (!props.vus.clinvarErrorMsg ||
                  props.vus.clinvarErrorMsg?.length === 0)
                  ? styles.disabled
                  : props.vus.rsid?.length > 0 && !props.vus.rsidDbsnpVerified
                  ? styles["unverified-rsid"]
                  : ""
              }`}
            >
              {props.vus.clinvarClassification?.length > 0 ? (
                <>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Classification:</div>
                    {props.vus.clinvarClassification}
                  </div>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Review status:</div>
                    {props.vus.clinvarClassificationReviewStatus}
                  </div>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Last evaluated:</div>
                    {props.vus.clinvarClassificationLastEval}
                  </div>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Canonical SPDI:</div>
                    {props.vus.clinvarCanonicalSpdi}
                  </div>
                </>
              ) : props.vus.clinvarErrorMsg?.length > 0 ? (
                <div className={styles.information}>
                  <div className={styles["info-title"]}>Error message:</div>
                  {props.vus.clinvarErrorMsg}
                </div>
              ) : (
                <div className={styles.information}>
                  No Clinvar entry found based on dbSNP's RSID
                </div>
              )}
            </div>
          </div>

          {/** DbSnp */}
          <div className={styles["info-container"]}>
            <div className={styles["external-ref-title-container"]}>
              <p className={styles["info-title"]}>dbSNP:</p>
              <Icon
                name="external-link"
                className={`${styles["external-link"]} ${styles.dbsnp} ${
                  props.vus.rsidDbsnpVerified
                    ? styles.active
                    : props.vus.rsid?.length > 0
                    ? styles["unverified-rsid"]
                    : styles.disabled
                }`}
                onClick={(e) => {
                  if (props.vus.rsidDbsnpVerified) {
                    e.stopPropagation();
                    openInNewWindow(
                      `https://www.ncbi.nlm.nih.gov/snp/${props.vus.rsid}`
                    );
                  }
                }}
              />
            </div>

            <div
              className={`${styles.info} ${
                props.vus.rsidDbsnpVerified
                  ? ""
                  : props.vus.rsid?.length > 0
                  ? styles["unverified-rsid"]
                  : styles.disabled
              }`}
            >
              {props.vus.rsid?.length > 0 ? (
                <>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>
                      Is RSID verified:
                    </div>
                    {props.vus.rsidDbsnpVerified.toString()}
                  </div>
                  {props.vus.rsidDbsnpVerified ? (
                    <div className={styles.information}>
                      <div className={styles["info-title"]}>RSID:</div>
                      {props.vus.rsid}
                    </div>
                  ) : (
                    <>
                      {props.vus.rsid !== "NORSID" && (
                        <div className={styles.information}>
                          <div className={styles["info-title"]}>
                            Suggested RSID:
                          </div>
                          <a
                            href={`https://www.ncbi.nlm.nih.gov/snp/${props.vus.rsid}`}
                            target="_blank"
                          >
                            {props.vus.rsid}
                          </a>
                        </div>
                      )}
                      <div className={styles.information}>
                        <div className={styles["info-title"]}>
                          Error message:
                        </div>
                        {props.vus.rsid === "NORSID"
                          ? "No RSID found."
                          : props.vus.rsidDbsnpErrorMsgs}
                      </div>
                    </>
                  )}
                </>
              ) : (
                <div className={styles.information}>
                  No valid RSID found for this variant.
                </div>
              )}
            </div>
          </div>
        </div>

        {/** ACMG rules */}
        <div className={styles["acmg-rules"]}>
          <p className={styles["info-title"]}>ACMG rules:</p>
          <div
            className={`${styles["acmg-rules-info"]} ${
              isAcmgEditMenuOpen ? styles["acmg-edit-open"] : ""
            }`}
          >
            <div className={styles["acmg-rules-edit"]}>
              <AcmgRulesEdit
                variantId={props.vus.id}
                variantAcmgRuleIds={props.vus.acmgRuleIds}
                allAcmgRules={props.acmgRules}
                vusService={props.vusService}
                onMenuAcmgRuleHover={(acmgRuleId?: number) =>
                  setAcmgRuleHover(acmgRuleId)
                }
                onEditIconClick={(isOpen) => setIsAcmgEditMenuOpen(isOpen)}
              />
            </div>
            {isAcmgEditMenuOpen && (
              <div className={styles["acmg-info"]}>
                <AcmgRuleInfo
                  acmgRule={props.acmgRules.find((r) => r.id === acmgRuleHover)}
                />
              </div>
            )}
          </div>
        </div>

        {/** Samples */}
        <div className={styles["samples-container"]}>
          <p className={styles["info-title"]}>Samples with this variant:</p>
          <div className={styles.samples}>
            {props.vus.samples.map((s) => (
              <div className={styles.sample}>
                <div className={styles.bullet}>{"\u25CF"}</div>
                <div>
                  <p onClick={() => openInNewWindow(`/sample/${s}`)}>{s}</p>{" "}
                </div>
              </div>
            ))}
          </div>
        </div>

        {/** Phenotypes */}
        <div className={styles["phenotypes-container"]}>
          <p className={styles["info-title"]}>
            Phenotypes of the above samples:
          </p>
          <p className={styles["info-description"]}>
            Checkout if there are any publications for this variant in relation
            the a phenotype by clicking on the book icon next to the phenotype.
          </p>
          <div className={styles.phenotypes}>
            {props.vus.phenotypes.map((p) => (
              <div className={styles["phenotype-container"]}>
                <div className={styles.bullet}>{"\u25CF"}</div>
                <div
                  className={styles.phenotype}
                  onClick={() =>
                    openInNewWindow(
                      `https://hpo.jax.org/app/browse/term/${p.ontologyId}`
                    )
                  }
                >
                  <b>{p.ontologyId}</b>: {p.name}
                </div>
                <Icon
                  name="publication"
                  className={styles["pub-icon"]}
                  onClick={() =>
                    openInNewWindow(
                      `/publication-phenotype-view/${props.vus.id}?rsid=${props.vus.rsid}&phenotype=${p.name}`
                    )
                  }
                />
              </div>
            ))}
          </div>
        </div>
      </div>
    </div>
  );
};

export default VusInfo;
