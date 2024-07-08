import { IAcmgRule } from "../../models/acmg-rule.model";
import { IClinvarUpdate } from "../../models/clinvar-updates.model";
import { IPhenotype } from "../../models/phenotype.model";
import { IPublicationPreview } from "../../models/publication-view.model";
import { ISampleToAddInfo } from "../../models/sample-to-add-info.model";
import { IUpdatedExternalRefData } from "../../models/updated-external-ref-data.model";
import { IVariantPublicationUpdates } from "../../models/variant-publication-updates";
import { INotVusSample, IVus, IVusSample } from "../../models/view-vus.model";
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

export interface ICheckFileForMultipleGenesRequest {
  vusFile: File;
}

export interface ICheckFileForMultipleGenesResponse {
  isSuccess: boolean;
  multipleGenes?: IVusGene[] | null;
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
  taskId: number;
}

export interface ICheckFileUploadStatusesRequest {
  taskIds: string[];
}

export interface IStatus {
  taskId: number;
  isSuccess?: boolean | null;
  filename?: string | null;
  areRsidsRetrieved?: boolean | null;
  isClinvarAccessed?: boolean | null;
  vusList?: IVUSSummary[] | null;
  noHpoTermPhenotypes?: INoHPOTermPhenotype[] | null;
}

export interface ICheckFileUploadStatusesResponse {
  statuses: IStatus[];
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

export interface IGetClinvarUpdatesRequest {
  clinvarId: number;
}

export interface IGetClinvarUpdatesResponse {
  isSuccess: boolean;
  clinvarUpdates: IClinvarUpdate[];
  datesWithUpdates?: string[] | null;
}

export interface IGetPublicationUpdatesRequest {
  variantId: number;
}

export interface IGetPublicationUpdatesResponse {
  isSuccess: boolean;
  variantPublicationUpdates: IVariantPublicationUpdates[];
  datesWithUpdates?: string[] | null;
}

export interface IAddPublicationsRequest {
  variantId: string;
  publicationUrls: string[];
}

export interface IAddPublicationsResponse {
  isSuccess: boolean;
  publications?: IPublicationPreview[];
}

export interface IDeleteVariantRequest {
  variantId: string;
}

export interface IDeleteVariantResponse {
  isSuccess: boolean;
}

export interface IAddSamplesRequest {
  variantId: number;
  samplesToAdd: ISampleToAddInfo[];
}

export interface IAddSamplesResponse {
  isSuccess: boolean;
  updatedSamples: IVusSample[];
  updatedNotVariantSamples: INotVusSample[];
  updatedPhenotypes: IPhenotype[];
}

export interface IAddNewSampleRequest {
  variantId: number;
  sampleToAdd: ISampleToAddInfo;
}

export interface IAddNewSampleResponse {
  isSuccess: boolean;
  updatedSamples?: IVusSample[];
  updatedNotVariantSamples?: INotVusSample[];
  updatedPhenotypes?: IPhenotype[];
}

export interface IRemoveSamplesRequest {
  variantId: number;
  sampleIdsToRemove: string[];
}

export interface IRemoveSamplesResponse {
  isSuccess: boolean;
  updatedSamples?: IVusSample[];
  updatedNotVariantSamples?: INotVusSample[];
  updatedPhenotypes?: IPhenotype[];
}

export interface IUpdateRsidRequest {
  variantId: number;
  newRsid: string;
}

export interface IUpdateRsidResponse {
  isSuccess: boolean;
  updatedExternalRefData?: IUpdatedExternalRefData | null
}