import { IHPOTerm } from "../services/vus/vus.dto";
import { IAcmgRule } from "./acmg-rule.model";
import { IVUSSummary } from "./publication-view.model";

export enum Genotype {
  Homozygous = "HOMOZYGOUS",
  Heterozygous = "HETEROZYGOUS",
}

export interface ISampleVariant {
  variantId: number;
  variant: IVUSSummary;
  genotype: Genotype;
  acmgRuleIds?: number[];
}

export interface ISample {
  sampleId: string;
  genomeVersion: string;
  fileUploadName: string;
  dateOfFileUpload: Date;
  variants: ISampleVariant[];
  phenotype?: IHPOTerm[] | null;
}

export interface ISampleSummary {
  sampleId: string;
  dateOfFileUpload: Date;
  numOfVariants: number;
}
