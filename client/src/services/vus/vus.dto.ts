import { IAcmgRule } from "../../models/acmg-rule.model";
import { IClinvarUpdate } from "../../models/clinvar-updates.model";
import { IVus } from "../../models/view-vus.model";
import {
  ISamplePhenotypeSelected,
  IVusGene,
  IVusGeneSelected,
} from "../../models/vus-file-upload.model";
import { IVUSSummary } from "../../models/vus-summary.model";
import { IVusUpload } from "../../models/vus-upload.model";

export interface IUploadVusRequest {
  vus: IVusUpload;
}

export interface IUploadVusResponse {
  isSuccess: boolean;
  vusList: IVus[];
}

export interface IStoreAndVerifyVusFileRequest {
  vusFile: File;
  multipleGenesSelection?: IVusGeneSelected[];
  samplePhenotypeSelection?: ISamplePhenotypeSelected[];
}

export interface INoHPOTermPhenotype {
  phenotype: string;
  samples: string[];
}

export interface IStoreAndVerifyVusFileResponse {
  isSuccess: boolean;
  areRsidsRetrieved: boolean;
  isClinvarAccessed: boolean;
  vusList: IVus[];
  multipleGenes?: IVusGene[] | null;
  noHpoTermPhenotypes: INoHPOTermPhenotype[];
}

export interface ILoadAllVusResponse {
  isSuccess: boolean;
  vusList?: IVUSSummary[] | null;
}

export interface IGetVusRequest {
  vusId: number;
}

export interface IGetVusResponse {
  isSuccess: boolean;
  vus: IVus;
  acmgRules: IAcmgRule[];
}

export interface IVerifyGeneRequest {
  geneName: string;
}

export interface IVerifyGeneResponse {
  isSuccess: boolean;
  geneId?: number | null;
}

export interface IGetAllAcmgRulesResponse {
  isSuccess: boolean;
  acmgRules: IAcmgRule[];
}

export interface IAddAcmgRuleRequest {
  variantId: number;
  ruleId: number;
}

export interface IAddAcmgRuleResponse {
  isSuccess: boolean;
}

export interface IRemoveAcmgRuleRequest {
  variantId: number;
  ruleId: number;
}

export interface IRemoveAcmgRuleResponse {
  isSuccess: boolean;
}

export interface IGetClinvarUpdatesRequest {
  clinvarId: number;
}

export interface IGetClinvarUpdatesResponse {
  isSuccess: boolean;
  clinvarUpdates: IClinvarUpdate[];
}
