import {
  IAddAcmgRuleRequest,
  IAddAcmgRuleResponse,
  IGetAllAcmgRulesResponse,
  IGetClinvarUpdatesRequest,
  IGetClinvarUpdatesResponse,
  IGetVusRequest,
  IGetVusResponse,
  ILoadAllVusResponse,
  IRemoveAcmgRuleRequest,
  IRemoveAcmgRuleResponse,
  IStoreAndVerifyVusFileRequest,
  IStoreAndVerifyVusFileResponse,
  IUploadVusRequest,
  IUploadVusResponse,
  IVerifyGeneRequest,
  IVerifyGeneResponse,
} from "./vus.dto";

export class VusService {
  async uploadVus(input: IUploadVusRequest): Promise<IUploadVusResponse> {
    let data = new FormData();

    const vusJsonData = input.vus ? JSON.stringify(input.vus) : "";

    // Append the JSON string as a blob to the FormData
    data.append("vus", vusJsonData);

    const result: IUploadVusResponse = await fetch(`/vus/upload`, {
      method: "POST",
      body: data,
      // headers: {
      // "Content-Type": "multipart/form-data",
      // },
      cache: "no-store", //TODO: is it needed?
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async storeAndVerifyVusFile(
    input: IStoreAndVerifyVusFileRequest
  ): Promise<IStoreAndVerifyVusFileResponse> {
    let data = new FormData();
    data.append("file", input.vusFile);

    const multipleGenesSelectionJsonData = input.multipleGenesSelection
      ? JSON.stringify(input.multipleGenesSelection)
      : "";

    // Append the JSON string as a blob to the FormData
    data.append("multipleGenesSelection", multipleGenesSelectionJsonData);

    const result: IStoreAndVerifyVusFileResponse = await fetch(`/vus/file`, {
      method: "POST",
      body: data,
      // headers: {
      // "Content-Type": "multipart/form-data",
      // },
      cache: "no-store", //TODO: is it needed?
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async loadAllVus(): Promise<ILoadAllVusResponse> {
    const result: ILoadAllVusResponse = await fetch(`/vus/view`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getVus(input: IGetVusRequest): Promise<IGetVusResponse> {
    const result: IGetVusResponse = await fetch(`/vus/view/${input.vusId}`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async verifyGene(input: IVerifyGeneRequest): Promise<IVerifyGeneResponse> {
    const result: IVerifyGeneResponse = await fetch(
      `/vus/gene/${input.geneName}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getAllAcmgRules(): Promise<IGetAllAcmgRulesResponse> {
    const result: IGetAllAcmgRulesResponse = await fetch(
      `/vus/all-acmg-rules`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async addAcmgRule(input: IAddAcmgRuleRequest): Promise<IAddAcmgRuleResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("variantId", input.variantId.toString());
    data.append("ruleId", input.ruleId.toString());

    const result: IAddAcmgRuleResponse = await fetch(`/vus/add-acmg-rule`, {
      method: "POST",
      body: data,
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async removeAcmgRule(
    input: IRemoveAcmgRuleRequest
  ): Promise<IRemoveAcmgRuleResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("variantId", input.variantId.toString());
    data.append("ruleId", input.ruleId.toString());

    const result: IRemoveAcmgRuleResponse = await fetch(
      `/vus/remove-acmg-rule`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getClinvarUpdates(
    input: IGetClinvarUpdatesRequest
  ): Promise<IGetClinvarUpdatesResponse> {
    const result: IGetClinvarUpdatesResponse = await fetch(
      `/vus/get_clinvar_updates/${input.clinvarId}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
        cache: "no-cache",
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const vusService = new VusService();
